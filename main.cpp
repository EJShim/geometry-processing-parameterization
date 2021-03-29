#include "tutte.h"
#include "lscm.h"
// #include <igl/read_triangle_mesh.h>
#include <igl/per_vertex_normals.h>
#include <igl/opengl/glfw/Viewer.h>
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
	mapper->SetScalarRange(1, 17);

	vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
	lut->SetUseBelowRangeColor(true);
	lut->SetBelowRangeColor(1, 1, 1, 1);
	mapper->SetLookupTable(lut);
	mapper->SetInputData(polydata);


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


int main(int argc, char *argv[])
{
	
	vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	iren->SetInteractorStyle( vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New() );
	vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
	iren->SetRenderWindow(renWin);
	vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
	renWin->AddRenderer(ren);



	//Read Mesh using VTK
	// vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
	// reader->SetFileName( "../data/animal.obj" );
	vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
	reader->SetFileName( "//192.168.0.113/Imagoworks/Data/confident/Mesh/IntraoralScan/DAEYOU-cut/train/2930/136.vtp" );
	reader->Update();
	std::cout << "OBJ File Reading Done! " << std::endl;



	vtkSmartPointer<vtkPolyData> polydata = reader->GetOutput();
	vtkSmartPointer<vtkActor> originalActor = MakeActor(polydata);
	// ren->AddActor(originalActor);





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
	
	



	Eigen::MatrixXd U_lscm, U_tutte, U;
	tutte(V,F,U_tutte);
	lscm(V,F,U_lscm);

	normalize(V);
	normalize(U_tutte);
	normalize(U_lscm);

	Eigen::MatrixXd V_tutte = (Eigen::MatrixXd(V.rows(),3) << U_tutte.col(0),Eigen::VectorXd::Zero(V.rows()),U_tutte.col(1)).finished();
	Eigen::MatrixXd V_lscm = (Eigen::MatrixXd(V.rows(),3) << U_lscm.col(0),Eigen::VectorXd::Zero(V.rows()),U_lscm.col(1)).finished();
	


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


	// auto lscmPoly = UpdateV(V_lscm, polydata);
	// vtkSmartPointer<vtkActor> lscmActor = MakeActor(lscmPoly);
	// lscmActor->SetPosition(3, 0, 0);
	// lscmActor->GetProperty()->SetColor(1, 0, 0);
	// ren->AddActor(lscmActor);



	// igl::opengl::glfw::Viewer viewer;
	// bool plot_parameterization = false;
	// const auto & update = [&]()
	// {
	// 	if(plot_parameterization)
	// 	{
	// 		// Viewer wants 3D coordinates, so pad UVs with column of zeros
	// 		viewer.data().set_vertices( V_tutte );
	// 		std::cout << "Plot Parameterization" << std::endl;
	// 	}else
	// 	{
	// 		viewer.data().set_vertices(V);
	// 	}


	// 	viewer.data().compute_normals();
	// 	viewer.data().set_uv(U*10);
	// };
	// viewer.callback_key_pressed = 
	// 	[&](igl::opengl::glfw::Viewer &, unsigned int key, int)
	// 	{
	// 		switch(key)
	// 		{
	// 		case ' ':
	// 			plot_parameterization ^= 1;
	// 			break;
	// 		case 'l':
	// 			U = U_lscm;
	// 			break;
	// 		case 't':
	// 			U = U_tutte;
	// 			break;
	// 		case 'C':
	// 		case 'c':
	// 			viewer.data().show_texture ^= 1;
	// 			break;
	// 		default:
	// 			return false;
	// 		}
	// 		update();
	// 		return true;
	// 	};

	// U = U_lscm;
	// viewer.data().set_mesh(V,F);
	// Eigen::MatrixXd N;
	// igl::per_vertex_normals(V,F,N);
	// viewer.data().set_colors(N.array()*0.5+0.5);
	// update();
	// viewer.data().show_texture = true;
	// viewer.data().show_lines = false;
	// viewer.launch();

	ren->ResetCamera();
	renWin->Render();
	iren->Start();



	return EXIT_SUCCESS;
}
