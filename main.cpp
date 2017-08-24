//
//  main.cpp
//  Random Forest
//
//  Created by zhang on 2017/7/16.
//  Copyright © 2017年 1. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include "stdlib.h"
#include "RandomForest.hpp"



int main(){
    //    int a,b, sum;
    //    cout << "a:" << endl;
    //    cin >> a;
    //    cout << "b:";
    //   cin >> b;
    //    sum = a + b;
    //    cout << "sum:" << a << "+" << b << "=" << sum << endl;
    //   return 0;
    
    // read the data from txt file
    /*
     ifstream file;
     char filename[512];
     cout << "File Name please:" <<endl;
     cin >> filename;
     file.open(filename,ios::in);
     if(file.fail()){
     cout<< "No file existing" <<endl;
     file.close();
     cin.get();
     cin.get();
     }
     else{
     
     }*/
    //int Bsize = 2000;
    
    //data: orginial sample
    
    //ofstream outfile;
    //outfile.open ("result.txt");
    std::ifstream file;
    file.open("/Users/zhang/Desktop/sim/sim3.txt");
    
    double data[Size][dim] = {-1};
    int i,j,k,h,Mi,Mj,deep=0,T=0;
    double d[2*dim];
    for(i = 0; i < 2*dim;i++){
        if(i%2 == 0){d[i] = 0;}
        else{d[i] = 1;}
    }
    
    /////*****  generate data (from R)  *****/////
    //default_random_engine rd; //time(NULL)
    //normal_distribution<double> normal(0,1);
    
    for(i = 0; i < Size; i++){
        for(j = 0; j< dim; j++){
            file >> data[i][j]; //= /*normal(rd)*/ randomBeta(2,3)
            std::cout << data[i][j] << std::endl;
            //outfile << data[i][j] << ", ";
        }
        //outfile << "." << endl;
    }
    //outfile.close();
    file.close();
    
    // notice the random forest struct
    //TNode* RTree;
    //RTree = (TNode*)malloc(sizeof(TNode)*1);
    TNode* *TList = new TNode*[N1*N2];
    
    double sample[Bsize][dim];
    
    //bootstrap
    //default_random_engine random;
    //uniform_int_distribution<int> unif1(0, 20000);
    //uniform_int_distribution<int> unif2(0, 10);
    //loop for the bootstrap
    //srand(time(NULL));
    for(Mi = 0;Mi < N1; Mi++){
        //generate the bootstrap sample
        for(k = 0; k < Bsize; k++){
            h = rand()%Size;
            for(j = 0; j < dim; j++){
                sample[k][j] = data[h][j];
            }
        }
        
        for(Mj = 0; Mj < N2; Mj++){
            TList[T] = new TNode;
            
            RandomForest(TList[T],sample,Bsize,d, 1, deep);
            T++;
        }
    }
    
    
    //obtain density
    
    int Cc = 0;
    double X[10000],Y[10000];
    
    std::ofstream outfile;
    outfile.open ("result.txt");
    
    for(i = 0; i < 10000; i++){
        X[i] = 1.0*i/10000;
        
        Y[i] = Marginal_Estiamte(X[i], TList, 0);
        outfile << X[i] << "," << Y[i] << std::endl;
        
        if(i>0 && Y[i] != Y[i-1]){
            Cc++;
        }
    }
    outfile.close();
    
    
    //std::cout << Marginal_Estiamte(0.5, TList, 5) << std::endl;
    
    
    std::cout << " H-error: " << Hellinger_Error(TList) << std::endl;
    
    Marginal_error(X,Y);
    
    std::cout << "marginal bin: " << Cc << std::endl;
    Eliminate_Forest(TList);
    return 0;
    
    
}
