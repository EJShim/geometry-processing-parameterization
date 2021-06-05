#include "tutte.h"
#include "lscm.h"
#include <igl/per_vertex_normals.h>
#include <Eigen/Core>
#include <string>
#include <iostream>
#include <vtkSmartPointer.h>
#include <vtkOBJReader.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkDataArray.h>
#include <vtkCell.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkOpenGLSphereMapper.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkAutoInit.h>
#include <vtkScalarsToColors.h>
#include <vtkLookupTable.h>
#include <vtkCamera.h>
#include <vtkPolyDataNormals.h>
#include <vtkUnsignedCharArray.h>
#include <vtkCurvatures.h>
#include <vtkSTLReader.h>
#include <vtkIdFilter.h>
#include <vtkFeatureEdges.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkStripper.h>
#include <vtkExtractEdges.h>


VTK_MODULE_INIT(vtkRenderingOpenGL2)
VTK_MODULE_INIT(vtkInteractionStyle)

void normalize(  Eigen::MatrixXd &U ){
	U.rowwise() -= U.colwise().mean().eval();
	U.array() /= 
	(U.colwise().maxCoeff() - U.colwise().minCoeff()).maxCoeff()/2.0;
}

vtkSmartPointer<vtkActor> MakeActor( vtkSmartPointer<vtkPolyData> polydata ){

	// vtkSmartPointer<vtkOpenGLSphereMapper> mapper = vtkSmartPointer<vtkOpenGLSphereMapper>::New();
	// mapper->SetRadius(0.01);
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	// mapper->SetScalarRange(1, 17);

	vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
	lut->SetUseBelowRangeColor(true);
	lut->SetBelowRangeColor(1, 1, 1, 1);
	lut->SetRange(1, 17);
	// mapper->SetLookupTable(lut);
	mapper->SetInputData(polydata);
	mapper->SetScalarRange(-0.5, 0.5);

	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);


	return actor;
}


vtkSmartPointer<vtkPolyData> UpdateV(Eigen::MatrixXd &V, vtkSmartPointer<vtkPolyData> polydata){
	vtkSmartPointer<vtkPolyData> result = vtkSmartPointer<vtkPolyData>::New();
	result->DeepCopy(polydata);

	if(result->GetPointData()->GetNormals()) result->GetPointData()->RemoveArray("Normals");
	// result->SetPoints(polydata->GetPoints());
	// result->SetPolys(polydata->GetPolys());

	for(vtkIdType i=0  ; i<result->GetNumberOfPoints() ; i++){
		double point[3] = { V(i,0), V(i,1), V(i,2)  };
		
		result->GetPoints()->SetPoint(i, point);
	}

	// if(result->GetPointData()->GetNormals() ) result->GetPointData()->RemoveArray("Normals");

	result->GetPoints()->Modified();

	return result;
}

Eigen::VectorXi GenerateBoundary(vtkSmartPointer<vtkPolyData> polydata){

	vtkSmartPointer<vtkIdFilter> idFilter = vtkSmartPointer<vtkIdFilter>::New();
	idFilter->SetInputData(polydata);
	idFilter->SetPointIdsArrayName("ids");
	idFilter->SetPointIds(true);
	idFilter->SetCellIds(false);

	vtkSmartPointer<vtkFeatureEdges> featureEdges = vtkSmartPointer<vtkFeatureEdges>::New();
	featureEdges->SetInputConnection(idFilter->GetOutputPort());
	featureEdges->BoundaryEdgesOn();
	featureEdges->FeatureEdgesOff();
	featureEdges->ManifoldEdgesOff();
	featureEdges->NonManifoldEdgesOff();

	vtkSmartPointer<vtkPolyDataConnectivityFilter> conFilter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
	conFilter->SetInputConnection(featureEdges->GetOutputPort());
	conFilter->SetExtractionModeToLargestRegion();


	vtkSmartPointer<vtkStripper> stripper = vtkSmartPointer<vtkStripper>::New();
	stripper->SetInputConnection(conFilter->GetOutputPort());
	stripper->JoinContiguousSegmentsOn();

	vtkSmartPointer<vtkCleanPolyData> cleanPoly = vtkSmartPointer<vtkCleanPolyData>::New();
	cleanPoly->SetInputConnection(stripper->GetOutputPort());
	
	vtkSmartPointer<vtkStripper> stripper2 = vtkSmartPointer<vtkStripper>::New();
	stripper2->SetInputConnection(cleanPoly->GetOutputPort());
	stripper2->JoinContiguousSegmentsOn();

	vtkSmartPointer<vtkCleanPolyData> cleanPoly2 = vtkSmartPointer<vtkCleanPolyData>::New();
	cleanPoly2->SetInputConnection(stripper2->GetOutputPort());
	cleanPoly2->Update();
	
	vtkSmartPointer<vtkPolyData> boundaryPoly = cleanPoly2->GetOutput();


	Eigen::VectorXi boundary( boundaryPoly->GetNumberOfPoints());
	
	for(int i=0 ; i<boundaryPoly->GetNumberOfPoints() ; i++){
		
		int tuple = boundaryPoly->GetPointData()->GetArray("ids")->GetTuple1(i);
		boundary(i) = tuple;		
	}	

	return boundary;
}


int main(int argc, char *argv[])
{

	if(argc < 2){
		argv[1] = "../Mandibular_orig.stl";
	}

	std::cout << argv[1] << std::endl;

	
	vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	iren->SetInteractorStyle( vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New() );
	vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
	iren->SetRenderWindow(renWin);
	vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
	ren->GetActiveCamera()->SetPosition(0, -1, 0);
	renWin->AddRenderer(ren);
	renWin->SetSize(512, 512);



	//Read Mesh using VTK
	// vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
	// reader->SetFileName( "../data/animal.obj" );
	vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
	reader->SetFileName( argv[1] );
	reader->Update();
	std::cout << "File Reading Done! " << std::endl;



	vtkSmartPointer<vtkPolyData> polydata = reader->GetOutput();	
	
	//Curvature

	vtkSmartPointer<vtkCurvatures> curvaturesFilter = vtkSmartPointer<vtkCurvatures>::New();
	curvaturesFilter->SetInputData(polydata);
	curvaturesFilter->SetCurvatureTypeToMinimum();
	// curvaturesFilter->SetCurvatureTypeToMaximum();
	// curvaturesFilter->SetCurvatureTypeToGaussian();
	// curvaturesFilter->SetCurvatureTypeToMean();
	curvaturesFilter->Update();


	polydata = curvaturesFilter->GetOutput();

	

	//Convert to Eigen Matrix
	Eigen::MatrixXd V( polydata->GetNumberOfPoints(), 3 );
	Eigen::MatrixXi F( polydata->GetNumberOfCells(), 3 );

	for(int i=0 ; i<polydata->GetNumberOfPoints() ; i++){
		double* point = polydata->GetPoint(i);

		// std::cout << point[0] << "," << point[1] << "," << point[2] << " // " << V.row(i) << std::endl;
		V(i, 0) = point[0];
		V(i, 1) = point[1];
		V(i, 2) = point[2];
	}

	for(int i=0 ; i<polydata->GetNumberOfCells() ; i++){
		vtkIdList* ids =  polydata->GetCell(i)->GetPointIds();		
		
		if(ids->GetNumberOfIds() != 3) std::cout << F.row(i);
		F(i, 0) = ids->GetId(0);
		F(i, 1) = ids->GetId(1);
		F(i, 2) = ids->GetId(2); 
	}

	//Generate edges using VTK
	vtkSmartPointer<vtkExtractEdges> edgeExtractor = vtkSmartPointer<vtkExtractEdges>::New();
	edgeExtractor->SetInputData(polydata);
	edgeExtractor->Update();
	vtkSmartPointer<vtkPolyData> edgePoly = edgeExtractor->GetOutput();
	Eigen::MatrixXi E( edgePoly->GetNumberOfCells(), 2 );
	for(int i=0 ; i<edgePoly->GetNumberOfCells() ; i++){
		vtkIdList* ids =  edgePoly->GetCell(i)->GetPointIds();				
		E(i, 0) = ids->GetId(0);
		E(i, 1) = ids->GetId(1);		
	}


	//Generate Boundary using VTK
	Eigen::VectorXi boundary =  GenerateBoundary(polydata);
	

	//Calculate Parameterization
	Eigen::MatrixXd U_tutte = tutte(V, F, E, boundary);

	//Make it 3D with normalization
	normalize(V);
	normalize(U_tutte);
	Eigen::MatrixXd V_tutte = (Eigen::MatrixXd(V.rows(),3) << U_tutte.col(0),Eigen::VectorXd::Zero(V.rows()),U_tutte.col(1)).finished();		
	

	//Render
	auto normalizedPoly =  UpdateV( V, polydata );
	vtkSmartPointer<vtkActor> normalizedActor = MakeActor(normalizedPoly);
	// normalizedActor->SetPosition(0, 0, 0);
	ren->AddActor(normalizedActor);


	auto tuttePoly = UpdateV(V_tutte, polydata);
	vtkSmartPointer<vtkActor> tutteActor = MakeActor(tuttePoly);
	tutteActor->SetPosition(2, 0, 0);
	tutteActor->GetProperty()->SetRepresentationToWireframe();
	tutteActor->GetProperty()->SetColor(1, 0, 0);
	ren->AddActor(tutteActor);



	ren->ResetCamera();
	renWin->Render();
	iren->Start();



	return EXIT_SUCCESS;
}
