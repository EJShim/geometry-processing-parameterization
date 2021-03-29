#include "tutte.h"
#include "lscm.h"
#include <igl/read_triangle_mesh.h>
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


void normalize(  Eigen::MatrixXd &U ){
	U.rowwise() -= U.colwise().mean().eval();
	U.array() /= 
	(U.colwise().maxCoeff() - U.colwise().minCoeff()).maxCoeff()/2.0;
}



int main(int argc, char *argv[])
{
	Eigen::MatrixXd U_lscm, U_tutte, U;


	//Read Mesh using VTK
	vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
	reader->SetFileName( (argc>1?argv[1]:"../data/beetle.obj") );
	reader->Update();
	std::cout << "OBJ File Reading Done! " << std::endl;



	vtkSmartPointer<vtkPolyData> polydata = reader->GetOutput();


	//Convert to Eigen Matrix
	Eigen::MatrixXd V( polydata->GetNumberOfPoints(), 3 );
	Eigen::MatrixXi F( polydata->GetNumberOfCells(), 3 );

	for(int i=0 ; i<polydata->GetNumberOfPoints() ; i++){
		double* point;
		polydata->GetPoint(i, point);

		// std::cout << point[0] << "," << point[1] << "," << point[2] << " // " << V.row(i) << std::endl;
		V(i, 0) = point[0];
		V(i, 1) = point[1];
		V(i, 2) = point[2];
	}

	std::cout << F.rows() << "," << F.cols() << std::endl;
	std::cout << polydata->GetNumberOfCells() << std::endl;


	for(int i=0 ; i<polydata->GetNumberOfCells() ; i++){
		vtkIdList* ids =  polydata->GetCell(i)->GetPointIds();		
		
		if(ids->GetNumberOfIds() != 3) std::cout << F.row(i);
		F(i, 0) = ids->GetId(0);
		F(i, 1) = ids->GetId(1);
		F(i, 2) = ids->GetId(2); 

	}
	
	igl::opengl::glfw::Viewer viewer;

	tutte(V,F,U_tutte);
	lscm(V,F,U_lscm);

	normalize(V);
	normalize(U_tutte);
	normalize(U_lscm);

	bool plot_parameterization = false;
	const auto & update = [&]()
	{
		if(plot_parameterization)
		{
		// Viewer wants 3D coordinates, so pad UVs with column of zeros
		viewer.data().set_vertices(
			(Eigen::MatrixXd(V.rows(),3)<<
			U.col(0),Eigen::VectorXd::Zero(V.rows()),U.col(1)).finished());
		}else
		{
		viewer.data().set_vertices(V);
		}
		viewer.data().compute_normals();
		viewer.data().set_uv(U*10);
	};
	viewer.callback_key_pressed = 
		[&](igl::opengl::glfw::Viewer &, unsigned int key, int)
	{
		switch(key)
		{
		case ' ':
			plot_parameterization ^= 1;
			break;
		case 'l':
			U = U_lscm;
			break;
		case 't':
			U = U_tutte;
			break;
		case 'C':
		case 'c':
			viewer.data().show_texture ^= 1;
			break;
		default:
			return false;
		}
		update();
		return true;
	};

	U = U_tutte;
	viewer.data().set_mesh(V,F);
	Eigen::MatrixXd N;
	igl::per_vertex_normals(V,F,N);
	viewer.data().set_colors(N.array()*0.5+0.5);
	update();
	viewer.data().show_texture = true;
	viewer.data().show_lines = false;
	viewer.launch();

	return EXIT_SUCCESS;
}
