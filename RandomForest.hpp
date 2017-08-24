//
//  RandomForest.hpp
//  Random Forest
//
//  Created by zhang on 2017/7/16.
//  Copyright © 2017年 1. All rights reserved.
//

#ifndef RandomForest_hpp
#define RandomForest_hpp

#include <stdio.h>
#define dim 5
#define Size 10000
#define E_size 1000000
#define theta 0.001
#define Epsilon 0.001
#define N1 10
#define N2 1


const int Bsize = 50000;



//int row,column;
//int r, f;



typedef struct TNode
{
    double density;
    double Division_position;
    int Division_cooordinate;
    struct TNode* left;
    struct TNode* right;
    double volume;
    double bound[2*dim];
}TNode;

int RandomForest(TNode *RT, double dat[Bsize][dim], double n, double D[2*dim], double volume, int Deep);
double Star_Discrepancy( double dat[Bsize][dim],double a[dim], double b[dim], double Row, int P);
double One_Star_Discrepancy(double dat[Bsize][dim], int P, double a, double b, double Row);
double Imbalance(double dat[Bsize][dim],double a, double b, int d, int k, double volume, double Row, int ndivide);
int* Imb(double dat[Bsize][dim], double a[dim], double b[dim], double Row, int ndivide);
double One_Dimension_Discrepancy(double dat[Bsize][dim], double a[dim], double b[dim], double Row);
double randomBeta( double alpha, double beta);
double Density_Estimate(double dat[dim], TNode* Forest[N1*N2]);
double Marginal_Estiamte(double dat, TNode* Forest[N1*N2], int P);
double Travel_Node(double dat, TNode* D, double sum, int P);
double Eliminate_Forest(TNode* Forest[N1*N2]);
void Eliminate_Node(TNode* D);
double Hellinger_Error(TNode* Forest[N1*N2]);
double Marginal_error(double x[10000], double y[10000]);




#endif /* RandomForest_hpp */
