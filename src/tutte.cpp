#include "tutte.h"
#include "math.h"
#include <vector>
#include <igl/boundary_loop.h>
#include <igl/edges.h>
#include <fstream>
#include <igl/min_quad_with_fixed.h>
#include <igl/repdiag.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>



Eigen::MatrixXd tutte( Eigen::MatrixXd V, Eigen::MatrixXi F)
{


	//Calculate Boundary, make it some shape
	Eigen::VectorXi boundary;
	igl::boundary_loop(F,boundary);

	return tutte(V, F, boundary);

	// Eigen::MatrixXd boundary_uv(boundary.size(),2);
	// // igl::map_vertices_to_circle(V,boundary,boundary_uv);
	// double pi=3.14159265359;
	// for (int i=0; i<boundary.rows(); i++){
	// 	boundary_uv.row(i)=Eigen::RowVector2d(cos(pi*2/boundary.size()*i),sin(pi*2/boundary.size()*i));		
	// }

	// // igl::harmonic(V,F,boundary,boundary_uv,0,U);
	// // return;	


	// //Calculate Edge
	// Eigen::MatrixXi E;
	// igl::edges(F,E);
	

	// Eigen::VectorXd diag=Eigen::VectorXd::Zero(V.rows());	
	// typedef Eigen::Triplet<double> T;
	// std::vector<T> list;
	// double tp;
	// for (int i=0; i<E.rows(); i++){
	// 	tp=1.0/(V.row(E(i,0))-V.row(E(i,1))).norm();
	// 	list.push_back(T(E(i,0),E(i,1),tp));
	// 	list.push_back(T(E(i,1),E(i,0),tp));
	// 	diag(E(i,0))-=tp;
	// 	diag(E(i,1))-=tp;
	// }
	// for (int i=0; i<V.rows(); i++)
	// 	list.push_back(T(i,i,diag(i)));

	// // std::cout << list.size() << std::endl;

	// Eigen::SparseMatrix<double> A;
	// A.resize(V.rows(),V.rows());
	// A.setFromTriplets(list.begin(), list.end()); 
	// const Eigen::VectorXd B_flat = Eigen::VectorXd::Zero(V.rows());
	// Eigen::SparseMatrix<double> Aeq = Eigen::SparseMatrix<double>();
	// Eigen::VectorXd Beq = Eigen::VectorXd();
	
	// Eigen::MatrixXd U_tutte;
	// U_tutte.resize(V.rows(),2);

	
	// igl::min_quad_with_fixed(A, B_flat, boundary, boundary_uv, Aeq, Beq, false, U_tutte  );



	// return U_tutte;
}


Eigen::MatrixXd tutte( Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::VectorXi boundary)
{
	

	Eigen::MatrixXd boundary_uv(boundary.size(),2);
	// igl::map_vertices_to_circle(V,boundary,boundary_uv);
	double pi=3.14159265359;
	for (int i=0; i<boundary.rows(); i++){
		boundary_uv.row(i)=Eigen::RowVector2d(cos(pi*2/boundary.size()*i),sin(pi*2/boundary.size()*i));		
	}

	// igl::harmonic(V,F,boundary,boundary_uv,0,U);
	// return;	


	//Calculate Edge
	Eigen::MatrixXi E;
	igl::edges(F,E);
	

	Eigen::VectorXd diag=Eigen::VectorXd::Zero(V.rows());	
	typedef Eigen::Triplet<double> T;
	std::vector<T> list;
	double tp;
	for (int i=0; i<E.rows(); i++){
		tp=1.0/(V.row(E(i,0))-V.row(E(i,1))).norm();
		list.push_back(T(E(i,0),E(i,1),tp));
		list.push_back(T(E(i,1),E(i,0),tp));
		diag(E(i,0))-=tp;
		diag(E(i,1))-=tp;
	}
	for (int i=0; i<V.rows(); i++)
		list.push_back(T(i,i,diag(i)));

	// std::cout << list.size() << std::endl;

	Eigen::SparseMatrix<double> A;
	A.resize(V.rows(),V.rows());
	A.setFromTriplets(list.begin(), list.end()); 
	const Eigen::VectorXd B_flat = Eigen::VectorXd::Zero(V.rows());
	Eigen::SparseMatrix<double> Aeq = Eigen::SparseMatrix<double>();
	Eigen::VectorXd Beq = Eigen::VectorXd();
	
	Eigen::MatrixXd U_tutte;
	U_tutte.resize(V.rows(),2);

	
	igl::min_quad_with_fixed(A, B_flat, boundary, boundary_uv, Aeq, Beq, false, U_tutte  );



	return U_tutte;
}