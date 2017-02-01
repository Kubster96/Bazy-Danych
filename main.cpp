#include<iostream>
#include<cstdio>
#include <cmath>
#include <fstream>
#include <algorithm>
using namespace std;

struct square{
    int vertexes[4];
};


 struct point{
  double x;
  double y;
  double z;
 } ;


 double g (double x, double y){

    return pow(x*x , 1.0/3.0);

 }

 bool sortComperator (point lhs, point rhs)
{
    if(fabs(lhs.x - rhs.x)<1e-6) return lhs.y > rhs.y;
    return lhs.x > rhs.x;
}


void fill_squares(square squares[], int n_of_squares, int n){

    int rows = 2*n;
    int columns = 2*n;


    for (int i=0;i<n;i++){
        for(int j=0;j<columns;j++){
            int k = i*rows + j;
            int l = i*(rows + 1) + j;
            squares[k].vertexes[3] = l;
            squares[k].vertexes[2] = l + 1;
            squares[k].vertexes[0] = l + columns + 1;
            squares[k].vertexes[1] = l + columns + 2;
        }
    }

    for(int i=2*n*n;i<n_of_squares;i++){
        squares[i].vertexes[3] = squares[i-n].vertexes[0];
        squares[i].vertexes[2] = squares[i].vertexes[3] + 1;
        squares[i].vertexes[0] = squares[i].vertexes[3] + n + 1;
        squares[i].vertexes[1] = squares[i].vertexes[0] + 1;
    }

}

int fill_bottom_vertex(int bottom_vertexes[],int n) {

    int counter = 0;
    for(int i=2*n*n+n;i<=2*n*n+2*n;i++){
        bottom_vertexes[counter]=i;
        counter++;
    }
    for(int i=2*n*n+3*n+1;i<=3*n*n+3*n;i+=n+1){
        bottom_vertexes[counter]=i;
        counter++;
    }
    return counter;
}



 void fill_points_coordinates( point points[], int size,int n){

       for(int i=0;i<size;i++){
        points[i].x=0.0;
        points[i].y=0.0;
        points[i].z=0.0;
       }

       double square_side=1.0/n;
       int point_index;

       for(int w=0;w<n+1;w++){                          // górny prostok¹t
         for(int k=0;k<2*n+1;k++){
            point_index = w * (2*n+1) + k;
            points[point_index].x=-1.0 + k*square_side;
            points[point_index].y = 1.0 - w*square_side;

         }

       }

      for(int w=0;w<n; w++){
        for(int k=0;k<n+1; k++){

            point_index = w * ( n +1) + k + (2*n +1)*(n+1);
            points[point_index].x= k*square_side;
            points[point_index].y = -square_side- w*square_side;

        }
      }

 }


 void fill_right_matrix (point points[], double right_matrix[], int size,int n){

    for(int i=0;i<size;i++) right_matrix[i]=0;

    double half_of_square_side = 1.0 / (2 * n);
    for(int i=1;i<2*n;i++){  //
        right_matrix[i] = ( g( points[i].x -half_of_square_side , points[i].y) + g(points[i].x + half_of_square_side, points[i].y) ) * half_of_square_side;
   }

    for (int i = size-2 ; i> size-n-1; i--){//
        right_matrix[i] = ( g( points[i].x -half_of_square_side , points[i].y) + g(points[i].x + half_of_square_side, points[i].y) ) * half_of_square_side;

    }
    for(int i=4*n+1; i < (2*n+1)*(n+1); i+=2*n+1){//
        right_matrix[i] = ( g( points[i].x , points[i].y-half_of_square_side) + g(points[i].x , points[i].y + half_of_square_side) ) * half_of_square_side;
    }


    for(int i=(2*n+1)*(n+1)+ n;i<size-1;i+=n+1){//
        right_matrix[i] = ( g( points[i].x , points[i].y-half_of_square_side) + g(points[i].x , points[i].y + half_of_square_side) ) * half_of_square_side;

    }
    for(int i=2*n+1; i <(2*n+1)*n ; i+=2*n+1){
        right_matrix[i] = ( g( points[i].x , points[i].y-half_of_square_side) + g(points[i].x , points[i].y + half_of_square_side) ) * half_of_square_side;
    }

    right_matrix[0] = ( g( points[0].x , points[0].y-half_of_square_side) + g(points[0].x +half_of_square_side , points[0].y) ) * half_of_square_side;
    right_matrix[2*n] = ( g( points[2*n].x , points[2*n].y-half_of_square_side) + g(points[2*n].x -half_of_square_side , points[2*n].y) ) * half_of_square_side;
    right_matrix[size-1] = ( g( points[size-1].x , points[size-1].y+half_of_square_side) + g(points[size-1].x -half_of_square_side , points[size-1].y) ) * half_of_square_side;

 }


int main(){

    int n;
    cout<<"Podaj n:"<<endl;
    cin>>n;

    int n_of_squares = 3*n*n;
    int n_of_vertexes = 3*n*n + 4*n + 1;
    int *bottom_vertexes = new int[2*n + 1];

    double psi[4][4] = {
         2.0/3.0,-1.0/6.0,-1.0/3.0,-1.0/6.0,
        -1.0/6.0, 2.0/3.0,-1.0/6.0,-1.0/3.0,
        -1.0/3.0,-1.0/6.0, 2.0/3.0,-1.0/6.0,
        -1.0/6.0,-1.0/3.0,-1.0/6.0, 2.0/3.0
    };


    double **left_matrix = new double*[n_of_vertexes];
    double *right_matrix = new double[n_of_vertexes];
    point *point_matrix = new point[n_of_vertexes];
    square *squares = new square [n_of_squares];

    for(int i=0;i<n_of_vertexes;i++){
        left_matrix[i] = new double[n_of_vertexes];
    }

    for(int i=0;i<n_of_vertexes;i++){
        for(int j=0;j<n_of_vertexes;j++){
            left_matrix[i][j] = 0.0;
        }
    }

    fill_squares(squares,n_of_squares,n);

    for(int i=0;i<n_of_squares;i++){
        for(int j=0;j<4;j++){
            for(int k=0;k<4;k++){
                left_matrix[squares[i].vertexes[j]][squares[i].vertexes[k]]+=psi[j][k];
            }
        }
    }

    int counter = fill_bottom_vertex(bottom_vertexes,n);

    for(int i=0;i<counter;i++){
        for(int j=0;j<n_of_vertexes;j++){
            left_matrix[bottom_vertexes[i]][j] = 0;
        }
        left_matrix[bottom_vertexes[i]][bottom_vertexes[i]] = 1;
    }


    fill_points_coordinates(point_matrix, n_of_vertexes, n);
    fill_right_matrix(point_matrix, right_matrix, n_of_vertexes, n);


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int N = n_of_vertexes;

   // Gaussian elimination with partial pivoting

    for(int i = 0; i < N; i++){
        // find pivot row and swap
        int max = i;
        for (int j = i + 1; j < N; j++)

        if(fabs(left_matrix[j][i]) > fabs(left_matrix[max][i])) max = j;

       //swap

        for(int k=0;k<N;k++){

            double temp;
            temp = left_matrix[max][k];
            left_matrix[max][k] = left_matrix[i][k];
            left_matrix[i][k] = temp;

        }
       //swap

        double temp2;

        temp2 = right_matrix[max];

        right_matrix[max] = right_matrix[i];
        right_matrix[i] = temp2;
       // pivot within b


        for(int j = i + 1; j < N; j++){

            double m = left_matrix[j][i] / left_matrix[i][i];
            if(fabs(m)<1e-6) m=0;
            right_matrix[j] -= right_matrix[i] * left_matrix[j][i] / left_matrix[i][i];

        }

        // pivot within A
        for(int j = i + 1; j < N; j++){

            double m = left_matrix[j][i] / left_matrix[i][i];

            if(fabs(m)<1e-6) m=0;

            for(int k = i+1; k < N; k++){

                left_matrix[j][k] -= left_matrix[i][k] * m;

            }
            left_matrix[j][i] = 0.0;

        }
    }

    double **result = new double *[N];

    for( int i = 0; i < N; ++i ){

        result[i] = new double [1];
        for ( int j = 0; j < 1; ++j) { result[i][j]=0.0; }

    }

    for (int j = N - 1; j >= 0; j--)
    {
        double t = 0.0;
        for (int k = j + 1; k < N; k++) { t += left_matrix[j][k] * result[k][0]; }
        result[j][0] = (right_matrix[j] - t) / left_matrix[j][j];
    }

    for(int i = 0; i<N; i++) {delete [] left_matrix[i]; }
    delete [] left_matrix;
    delete [] right_matrix;
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    fstream plik( "dane.dat", ios::out | ios::trunc );
    plik.close();
    plik.open( "dane.dat", std::ios::in | std::ios::out );
    if( plik.good() == true ){
        for(int i =0 ; i < N ; i++){
             point_matrix[i].z = result[i][0];
         }
         delete [] result;

         sort(point_matrix,point_matrix+N,sortComperator);
         double flag = point_matrix[0].x;
         for (int i =0 ; i < N-1 ; i++)
         {
             if(fabs(flag - point_matrix[i+1].x) > 1e-6)
             {
                 flag = point_matrix[i+1].x;
                 plik << point_matrix[i].x << "  " << point_matrix[i].y << "  " << point_matrix[i].z << "\n\n";
             }
             else
             {
             plik << point_matrix[i].x << "  " << point_matrix[i].y << "  " << point_matrix[i].z << "\n";
             }
         }
        plik << point_matrix[N-1].x << "  " << point_matrix[N-1].y << "  " << point_matrix[N-1].z << "\n\n";
         plik.close();
         delete [] point_matrix;
     }
     
    return 0;

}
