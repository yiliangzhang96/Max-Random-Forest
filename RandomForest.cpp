//
//  RandomForest.cpp
//  Random Forest
//
//  Created by zhang on 2017/7/16.
//  Copyright © 2017年 1. All rights reserved.
//

#include "RandomForest.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <math.h>
#include "stdlib.h"

//#define STACKSIZE 419430400
//#define ndivide 3
#define sub 3
std::ofstream outfile2;

//int Row;
//int deepth = 0;

// original rule for spliting
// a: left bound   b: right bound   d: dimension  k:number of (quantile) volume: volume of the cube
double Imbalance(double dat[Bsize][dim],double a, double b, int d, int k, double volume, double Row, int ndivide){
    int i;
    double p=0,q = 0;
    for(i = 0; i < Row; i++){
        //for(j = 0; j < column; j++){
        //   Imb[d] = (dat[i][d]-a)/b;
        //}
        if(((dat[i][d]-a)/b) < a + k*(b-a)/ndivide){p++;}
        else{q++;}
    }
    //std::cout << p << " " << q << " " << volume << " " << (p*k/ndivide + q*(1-k/ndivide)) << " ";
    //return fabs((p*k/ndivide-q*(1-k/ndivide))/((p*k/ndivide + q*(1-k/ndivide))*volume));
    return fabs((p*ndivide/k-q*(ndivide/(ndivide - k)))/((p*ndivide/k + q*(ndivide/(ndivide - k)))*(pow(-log2(b-a),2)+1)));
    //return (pow(p,2)*ndivide/k + pow(q,2)*ndivide/(ndivide - k) - pow(p+q,2));
}

int* Imb(double dat[Bsize][dim], double a[dim], double b[dim], double Row, int ndivide, int S[sub]){
    int i, j, s, num;
    double n;
    int res[2];
    res[0] = -1;
    res[1] = -1;
    
    
    int P;
    double test[ndivide];
    double c[ndivide];
    double prod, max;
    
    num = 0;
    max = 0;
    
    //generate vector test
    //i loop for test[100];  s loop for dat[];  sj loop for element in dat[][]
    for(P = 0; P<sub; P++){
        for(i = 1; i < ndivide; i++){
            test[i] = 1.00*i/ndivide;
            num = 0;
            for(s = 0; s<Row; s++){
                
                //k = 0;
                
                if(((dat[s][S[P]]-a[S[P]])) < b[S[P]]*test[i]){
                    num++;
                }
            }
            
            prod = test[i];
            c[i] = fabs((num/Row)-prod);
            //std::cout << "ci: " << c[i] << std::endl;
            if(c[i]>max){max = c[i]; res[0] = S[P]; res[1] = i;}
        }
    }
    
    
    
    /*for(i = 0; i<dim; i++){
     n =One_Star_Discrepancy(dat, i, a[i], b[i], Row);
     //std::cout << "N: " << n << std::endl;
     if(n > max){
     max = n;
     res[0] = i;
     }
     }*/
    std::cout << " dim: " << res[0] << " posit: "<< res[1] << std::endl;
    return res;
    
}


double One_Star_Discrepancy(double dat[Bsize][dim], int P, double a, double b, double Row){
    double test[100];
    double c[100];
    int i, j, s, num;
    double prod, max;
    
    num = 0;
    max = 0;
    
    //generate vector test
    //i loop for test[100];  s loop for dat[];  sj loop for element in dat[][]
    for(i = 0; i<100; i++){
        test[i] = (i+1)*0.01;
        num = 0;
        for(s = 0; s<Row; s++){
            
            //k = 0;
            
            if(((dat[s][P]-a)) < b*test[i]){
                num++;
            }
        }
        
        prod = test[i];
        c[i] = fabs((num/Row)-prod);
        //std::cout << "ci: " << c[i] << std::endl;
        if(c[i]>max){max = c[i];}
    }
    
    
    return max;
}


double One_Dimension_Discrepancy(double dat[Bsize][dim], double a[dim], double b[dim], double Row){
    int i;
    double n, max = 0;
    for(i = 0; i<dim; i++){
        n =One_Star_Discrepancy(dat, i, a[i], b[i], Row);
        //std::cout << "N: " << n << std::endl;
        if(n > max){
            max = n;
        }
    }
    
    return max;
}


// a,b:boundary for the cube;  Row: number of rows needed in calculation
double Star_Discrepancy(double dat[Bsize][dim],double a[dim], double b[dim], double Row, int P){
    
    
    
    
    double test[50000][dim];
    double c[50000];
    int i, j, k, s, num;
    double prod, max;
    
    
    //std::cout << & num << std::endl;
    //num = 0;
    max = 0;
    
    //generate vector test
    //i loop for test[20000];  j loop to initialize test;  s loop for dat[];  sj loop for element in dat[][]
    for(i = 0; i<50000; i++){
        //srand(time(NULL));
        for(j = 0; j<dim; j++){
            //k = 0;
            //std::default_random_engine random(time(NULL));
            //std::uniform_real_distribution<double> distribution(0,1);
            test[i][j] = (rand()%999)+1;
            test[i][j] = test[i][j]/1000;
        }
        
        num = 0;
        for(s = 0; s<Row; s++){
            
            k = 0;
            
            for(j = 0; j<dim; j++){
                //std::default_random_engine random;
                //std::uniform_real_distribution<double> distribution(0,1);
                //test[i][j] = distribution(random);
                
                if(((dat[s][j]-a[j])/b[j]) < test[i][j]){
                    k++;
                }
            }
            if(k==dim){
                num++;
            }
        }
        prod = 1;
        for(j = 0; j< dim; j++){
            prod = prod * test[i][j];
        }
        c[i] = fabs((num/Row)-prod);
        if(c[i]>max){max = c[i];}
    }
    
    //std::cout << "MAX: " << max <<std::endl;
    return max;
    
}





//RT:tree  dat:data that going to generate branch  n:number of data inside   D:cube information
int RandomForest(TNode *RT, double dat[Bsize][dim], double n, double D[2*dim], double volume, int Deep)
{
    
    int i,j,k;
    double max = 0, temp;
    int Nm = -1, NL, NR;
    double Im = -1;
    int s[sub];
    
    int ndivide = 3;
    
    //if(n/(Bsize*volume) > -1){ndivide = 3;}
    //else{ndivide = 4;}
    
    if(n==0){
        RT->density = 0;
        std::cout << "n==0! Depth = " << Deep << "; " << "v = "<< volume << "; " << "Density = " << RT->density << std::endl;
        return 0;
    }
    
    
     /*if(n/(Bsize*volume) > 3*35){
     RT->density = n/(Bsize*volume);
     
     for(i = 0;i < 2*dim; i++){
     RT->bound[i] = D[i];
     }
     RT->volume = volume;
     std::cout << "n==1! Depth = " << Deep << "; " << "v = "<< volume << "; " << "Density = " << RT->density << std::endl;
     
     if(RT->bound[0] < 0.3 && 0.3 < RT->bound[1] ){
     
     }
     
     
     return 0;
     
     }*/
    
    
    Deep++;
    //parameter manipulation
    double A[dim], B[dim];
    for(i = 0; i < 2*dim; i++){
        if(i%2 == 0){
            A[i/2] = D[i];
        }
        else{
            B[(i-1)/2] = D[i];
        }
    }
    
    RT->density = -1;
    //subsample
    //std::default_random_engine random(time(NULL));
    //std::uniform_int_distribution<int> unif3(0,dim-1);
    for(i = 0; i< sub; i++){
        s[i] = rand()%dim;
    }
    
    
    
    std::cout << "    First of all, v = " << volume << std::endl;
    
    //if n = 0

    
    /*if(n==1){
        RT->density = 1/(Bsize*volume);
        
        for(i = 0;i < 2*dim; i++){
            RT->bound[i] = D[i];
        }
        RT->volume = volume;
        std::cout << "n==1! Depth = " << Deep << "; " << "v = "<< volume << "; " << "Density = " << RT->density << std::endl;
        return 0;
        
    }*/
    
    //find the imbalance information
    for(k = 0; k < sub; k++){
     
     for(i = 1; i < ndivide; i++){
     //std::cout << n << "haha";
     temp = Imbalance(dat, D[2*s[k]], D[2*s[k]+1], s[k], i, volume, n, ndivide);
     std::cout << temp << std::endl;
     
              if(max < temp && (D[2*s[k]+1] - D[2*s[k]] > 0*pow(0.5,4))){ //Imbalance(dat, D[2*s[k]-1], D[2*s[k]], s[k], i, volume, n)
                  max = temp;//Imbalance(dat, D[2*s[k]-1], D[2*s[k]], s[k], i, volume, n);
                 Nm = s[k]; //split dimension
                 Im = i; //split position
             }
         }
     }
    
    /*int* S = Imb(dat, A, B, n, ndivide,s);
    
    Nm = S[0];
    Im = S[1];
     */
    
    
    if(Nm == -1){
        
        RT->density = n/(Bsize*volume);
        
        for(i = 0;i < 2*dim; i++){
            RT->bound[i] = D[i];
        }
        RT->volume = volume;
        
        std::cout << "Even! Depth = " << Deep << "; " << "v = "<< volume << "; " << "Density = " << RT->density << std::endl;
        return 0;
    }
    
    double threshold = theta * pow(Bsize,0.5)/n;
    
    //One_Star_Discrepancy(dat, Nm, A[Nm], B[Nm], n);
    //Star_Discrepancy(dat, A, B, n, Nm);
    //std::cout << "  N: " << n <<"  One: " << One_Dimension_Discrepancy(dat, A, B, n) << "  Star: " << Star_Discrepancy(dat, A, B, n, Nm) << "  Threshold: " << threshold << std::endl;
    
    
    
    //////*****      Split Condition     *****///////(One_Star_Discrepancy(dat, Nm, A[Nm], B[Nm], n) > threshold )
    if((Deep < 20) && ((threshold < Epsilon) || (One_Dimension_Discrepancy(dat, A, B, n) > threshold ) || (Star_Discrepancy(dat, A, B, n, Nm) > threshold  ))){
        std::cout << "split!";
        
        
        double LPass[Bsize][dim];
        double RPass[Bsize][dim];
        
        // Va split absolute position
        double Va = (Im * D[2*Nm + 1] + (ndivide - Im) * D[2*Nm])/ndivide;
        NL = NR = 0;
        for(i = 0; i < n; i++){
            if(dat[i][Nm] < Va){
                for(j = 0; j < dim; j++){
                    LPass[NL][j] = dat[i][j];
                }
                NL++;
            }
            else{
                for(j = 0; j < dim; j++){
                    RPass[NR][j] = dat[i][j];
                }
                NR++;
            }
        }
        /*
        if(n < 5 && (NR == 0 || NL == 0)){
            RT->density = n/(Bsize*volume);
            //RT->bound[0] = D[2*Nm];
            for(i = 0;i < 2*dim; i++){
                RT->bound[i] = D[i];
            }
            RT->volume = volume;
            std::cout << "Leaf with n ->> 0! Depth = " << Deep << "; " << "Density = " << RT->density << std::endl;
            
            return 0;
        }*/
        
        
        //std::cout << "  One: " << One_Dimension_Discrepancy(dat, A, B, n) << "  Star: " << Star_Discrepancy(dat, A, B, n, Nm) << "  Threshold: " << threshold << std::endl;
        //build left child; right child
        TNode* left;
        TNode* right;
        
        left = (TNode*)malloc(sizeof(TNode)*1);
        right = (TNode*)malloc(sizeof(TNode)*1);
        
        left->left = NULL;
        left->right = NULL;
        right->left = NULL;
        right->right = NULL;
        
        RT->left = left;
        RT->right = right;
        
        //NL to calculate the sample in left child;   LPass to be passed into left

        
        /*for(i = 0; i < Bsize; i++){
         for(j=0;j<dim;j++){
         LPass[i][j] = -1;
         RPass[i][j] = -1;
         }
         }*/
        
        if(Deep >= 23){
            
        }
        
    
            
        //0: coordinate;    1: position
        RT->Division_cooordinate = Nm;
        RT->Division_position = Va;
        for(i = 0;i < 2*dim; i++){
            RT->bound[i] = D[i];
        }
        //RT->bound[0] = D[2*Nm];
        //RT->bound[1] = D[2*Nm + 1];
        //RT->marginal_weight = volume/(D[2*Nm+1]-D[2*Nm]);
        std::cout << " At " << Nm << std::endl;
        
        
        //left branch
        //std::cout << "left: " << volume << " Im:" << Im << " prod: " << Im/ndivide << " In: " << volume*Im/ndivide << std::endl;
        double LD[2*dim];
        for(i=0;i<2*dim;i++){LD[i] = D[i];}
        LD[2*Nm+1] = Va;

        
        //?
        RandomForest(left, LPass, NL, LD, volume*Im/ndivide, Deep);
        
        
        //right branch
        //std::cout << "right: " << volume << " Im:" << Im << " prod: " << (1-(Im/ndivide)) << " In: " << volume*(1-(Im/ndivide)) << std::endl;
        double RD[2*dim];
        for(i=0;i<2*dim;i++){RD[i] = D[i];}
        RD[2*Nm] = Va;
        RandomForest(right, RPass, NR, RD, volume*(1-(Im/ndivide)), Deep);
    }
    
    //output density
    else{

        RT->density = n/(Bsize*volume);
        //RT->bound[0] = D[2*Nm];
        for(i = 0;i < 2*dim; i++){
            RT->bound[i] = D[i];
        }
        RT->volume = volume;
        std::cout << "Leaf! Depth = " << Deep << "; " << "Density = " << RT->density << std::endl;
    }
    
    return 0;
}









double randomBeta( double alpha, double beta)
{
    /*Johnk's beta generator*/
    double u, v;
    double x, y;
    
    //std::uniform_real_distribution<double> rb(0,5);
    do
    {
        //std::default_random_engine random1;
        u = rand() %3000;
        std::cout << u/1000 << std::endl;
        //std::default_random_engine random2;
        v = rand() %3000;
        std::cout << v/1000 << std::endl;
        x=pow(u/1000,1/alpha);
        y=pow(v/1000,1/beta);
    } while (x+y>1);
    return x/(x+y);
}



double Density_Estimate(double dat[dim], TNode* Forest[N1*N2])
{
    int i;
    double density=0;
    for (i = 0;i < N1*N2; i++){
        TNode* T = Forest[i];
        while(T->left != NULL){
            if(dat[T->Division_cooordinate] < T->Division_position){
                T = T->left;
            }
            else{
                T = T->right;
            }
        }
        
        density = density + T->density;
        
    }
    
    
    return density/(N1*N2);
}

double Marginal_Estiamte(double dat, TNode* Forest[N1*N2], int P){
    int i;
    double M = 0, marg[N1*N2];
    for(i = 0; i < N1*N2; i++){
        TNode* T = Forest[i];
        //sum = 0
        marg[i] = Travel_Node(dat, T, 0, P);
        //if(i == 9){
        M = M + marg[i];
        //}
        if(dat == 0.001 || dat == 0.998){
            
        }
    }
    
    return M/(N1*N2);
}

double Eliminate_Forest(TNode* Forest[N1*N2]){
    int i;
    outfile2.open ("calculate.txt");
    for(i = 0; i < N1*N2; i++){
        TNode* T = Forest[i];
        Eliminate_Node(T);
    }
    outfile2.close();
    return 0;
}

void Eliminate_Node(TNode* D){
    if(D->left->density != -1){
        free(D->left);
        
        
        outfile2 << 1 << std::endl;
        
    }
    else{
        Eliminate_Node(D->left);
    }
    
    if(D->right->density != -1){
        free(D->right);
    }
    else{
        Eliminate_Node(D->right);
    }
    
    free(D);
}





double Travel_Node(double dat, TNode* D, double sum, int P){
    
    if(D->density == -1){
        
        if(D->Division_cooordinate == P){
            
            if(dat < (D->Division_position)){
                sum = Travel_Node(dat, D->left, sum, P);
            }
            
            else{
                sum = Travel_Node(dat, D->right, sum, P);
            }
        }
        
        else{
            sum = Travel_Node(dat, D->left, sum, P);
            sum = Travel_Node(dat, D->right, sum, P);
        }
    }
    
    else{
        
        if(D->density == 0){
            
            sum = sum;
        }
        else{
            
            
            sum = sum + (D->density)*(D->volume)/(D->bound[2*P+1] - D->bound[2*P]);
            
        }
    }
    
    return sum;
}

double Star_Discrepancy2(double dat[Bsize][dim],double a[dim], double b[dim], double Row, int P){
    
    double test[100];
    double c[dim][100];
    int i, j, k, s, num;
    int prod, max;
    
    num = 0;
    max = 0;
    
    //generate vector test
    //i loop for test[100];  s loop for dat[];  sj loop for element in dat[][]
    for(i = 0; i<100; i++){
        test[i] = (i+1)/100;
        for(k = 0; k<dim; k++){
            for(s = 0; s<Row; s++){
                
                //k = 0;
                
                if(((dat[s][P]-a[P])/b[P]) < test[i]){
                    num++;
                }
            }
            
            prod = test[i];
            c[k][i] = fabs((num/Row)-prod);
            if(c[k][i]>max){max = c[k][i];}
        }
        
    }
    
    
    return max;
    
    
}

double Marginal_error(double x[10000], double y[10000]){
    int i,j;
    double f, E = 0, F = 0;
    for(i = 0; i < 10000;i++){
        f = (cos(x[i]*3.14*3.5)+1)/(sin(3.14*3.5)/(3.14*3.5) + 1);
        E = E + pow(f*y[i],0.5);
        F = F + pow(f-y[i],2);
    }
    std::cout << "Marginal_HD: " << 1 - E/10000 << std::endl;
    std::cout << "Marginal_MISE: " << F/10000 << std::endl;
    
    return 0;
}



double Hellinger_Error(TNode* Forest[N1*N2]){
    int i,j;
    double dat[dim];
    double E = 0, F = 0;
    double f,g;
    
    for(i = 0; i < E_size; i++){
        f = 1;
        for(j = 0; j < dim; j++){
            dat[j] = (rand()%10000)*0.0001;
            f = f*(cos(dat[j]*3.14*3.5)+1)/(sin(3.14*3.5)/(3.14*3.5) + 1);
            //f = f*(dat[j]*(1-dat[j])*(1-dat[j])*12);
        }
        
        //std::cout << " f: " << f << " Est:" << Density_Estimate(dat, Forest) << std::endl;
        g = Density_Estimate(dat, Forest);
        E = E + pow(f*g,0.5);
        F = F + pow(f - g,2);
    }
    
    std::cout << "ADG: " << F/E_size << std::endl;
    return 1 - E/E_size;
}
