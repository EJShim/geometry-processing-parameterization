#include "tutte.h"
#include "math.h"
#include <vector>
#include <igl/boundary_loop.h>
#include <igl/edges.h>
#include <fstream>
#include <igl/min_quad_with_fixed.h>
#include <igl/repdiag.h>

void tutte( const Eigen::MatrixXd & V, const Eigen::MatrixXi & F, Eigen::MatrixXd & U)
{
	std::vector<int> L;
	igl::boundary_loop(F,L);

	std::cout <<  L.size() << std::endl;
	Eigen::MatrixXi E;
	igl::edges(F,E);


	double pi=3.14159265359;
	Eigen::VectorXi b(L.size());
	Eigen::MatrixXd bc(L.size(),2);
	for (int i=0; i<L.size(); i++){
		bc.row(i)=Eigen::RowVector2d(cos(pi*2/L.size()*i),sin(pi*2/L.size()*i));
		b(i)=L[i];
	}

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

	std::cout << list.size() << std::endl;

	Eigen::SparseMatrix<double> A;
	A.resize(V.rows(),V.rows());
	A.setFromTriplets(list.begin(), list.end()); 
	const Eigen::VectorXd B_flat = Eigen::VectorXd::Zero(V.rows());
	Eigen::SparseMatrix<double> Aeq = Eigen::SparseMatrix<double>();
	Eigen::VectorXd Beq = Eigen::VectorXd();

	U.resize(V.rows(),2);

	std::cout << "A : " << A.rows() << "," << A.cols() <<  std::endl;
	std::cout << "B_flat : " << B_flat.rows() << "," << B_flat.cols() <<  std::endl;
	std::cout << "b : " << b.rows() << "," << b.cols() <<  std::endl;
	std::cout << "bc : " << bc.rows() << "," << bc.cols() <<  std::endl;
	std::cout << "Aeq : " << Aeq.rows() << "," << Aeq.cols() <<  std::endl;
	std::cout << "Beq : " << Beq.rows() << "," << Beq.cols() <<  std::endl;
	igl::min_quad_with_fixed(A, B_flat, b, bc, Aeq, Beq, false, U  );
}