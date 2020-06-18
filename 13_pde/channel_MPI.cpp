#include<iostream>
#include<vector>
#include<string>
#include<iterator>
#include<fstream>
#include<cmath>
#include"mpi.h"

using namespace std;

vector<vector<double> > build_up_b(double rho,double dt,double dx,double dy,vector<vector<double> > u,vector<vector<double> > v){
	vector< vector<double> > b(u.size(),vector<double>(u.at(0).size(),0));
	for(int i=1;i<b.size()-1;i++){
		for(int j=1;j<b.at(0).size()-1;j++){
			b[i][j] = (double)(rho*(1/dt*((u[i][j+1]-u[i][j-1]) / (2*dx) + (v[i+1][j]-v[i-1][j]) /(2*dy)) - pow(((u[i][j+1]-u[i][j-1]) / (2*dx)),2) - 2*((u[i+1][j]-u[i-1][j]) / (2*dy)*(v[i][j+1]-v[i][j-1])/(2*dx)) - pow(((v[i+1][j] - v[i-1][j]) / (2*dy)),2)));	
		}
	}
	int n = b.at(0).size();
	for(int i=1;i<b.size()-1;i++){
		b[i][n-1] = (double)(rho*(1/dt*((u[i][0]-u[i][n-2]) / (2*dx) + (v[i+1][n-1] - v[i-1][n-1]) / (2*dy)) - pow(((u[i][0]-u[i][n-2]) / (2*dx)),2) - 2*((u[i+1][n-1] - u[i-1][n-1]) / (2*dy) * (v[i][0]-v[i][n-2]) / (2*dx)) - pow(((v[i+1][n-1]-v[i-1][n-1]) / (2*dy)),2)));
	}

	for(int i=1;i<b.size()-1;i++){
		b[i][0] = (double)(rho*(1/dt*((u[i][1]-u[i][n-1]) / (2*dx) + (v[i+1][0] - v[i-1][0]) / (2*dy)) - pow(((u[i][1]-u[i][n-1]) / (2*dx)),2) - 2*((u[i+1][0]-u[i-1][0]) / (2*dy) * (v[i][1]-v[i][n-1]) / (2*dx)) - pow(((v[i+1][0]-v[i-1][0]) / (2*dy)),2)));
	}
	
	return b;

}


vector<vector<double> > pressure_poisson_periodic(vector<vector<double> > p,vector<vector<double> > b,double dx,double dy,int nit){
	vector< vector<double> > pn(p.size(),vector<double>(p.at(0).size(),0));
	
	for(int i=0;i<p.size();i++){
		for(int j=0;j<p.at(0).size();j++){
			pn[i][j] = p[i][j];
		}
	}

	for(int k=0;k<nit;k++){
	
		for(int i=1;i<p.size()-1;i++){
			for(int j=1;j<p.at(0).size()-1;j++){
				p[i][j] = (((pn[i][j+1]+pn[i][j-1])*pow(dy,2)+(pn[i+1][j]+pn[i-1][j])*pow(dx,2))/(2*(pow(dx,2)+pow(dy,2)))-pow(dx,2)*pow(dy,2)/(2*(pow(dx,2)+pow(dy,2)))*b[i][j]);
			}
		}
	
		int n = p.at(0).size();
		for(int i=1;i<p.size()-1;i++){
			p[i][n-1] = (((pn[i][0]+pn[i][n-2])*pow(dy,2)+(pn[i+1][n-1]+pn[i-1][n-1])*pow(dx,2))/(2*(pow(dx,2)+pow(dy,2)))-pow(dx,2)*pow(dy,2)/(2*(pow(dx,2)+pow(dy,2)))*b[i][n-1]);
		}
	
		for(int i=1;i<p.size()-1;i++){
			p[i][0] = (((pn[i][1]+pn[i][n-1])*pow(dy,2)+(pn[i+1][0]+pn[i-1][0])*pow(dx,2))/(2*(pow(dx,2)+pow(dy,2)))-pow(dx,2)*pow(dy,2)/(2*(pow(dx,2)+pow(dy,2)))*b[i][0]);
		}

		for(int j=0;j<p.at(0).size();j++){
			p[n-1][j] = p[n-2][j];
			p[0][j] = p[1][j];
		}		
	
	}
	return p;
}
int main(int argc,char** argv){
	MPI_Init(&argc,&argv);
	int size,rank;
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	int nx = 41;
	int ny = 41;
	int nt = 10;
	int nit = 50;
	int c = 1;
	double dx = (double) 2/(nx-1);
	double dy = (double) 2/(ny-1);
	vector<double> x(nx);
	vector<double> y(ny);
	double dif = 0;
	for (int i=0;i<nx;i++){
		x[i] = dif;
		dif += dx;
	}
	dif = 0;
	for (int i=0;i<ny;i++){
		y[i] = dif;
		dif += dy;
	}
	
	vector< vector<double> > X(nx,x);
	vector< vector<double> > Y(ny,vector<double>(ny));
	for(int i=0;i<ny;i++){
		for(int j=0;j<ny;j++){
			Y[i][j] = y[i];
		}
	}
	
	double rho = 1;
	double nu = 0.1;
	double F = 1;
	double dt = 0.01;
	vector<vector<double> > u(ny, vector<double>(nx, 0));
	vector<vector<double> > un(ny, vector<double>(nx, 0));
	vector<vector<double> > v(ny, vector<double>(nx, 0));
	vector<vector<double> > vn(ny, vector<double>(nx, 0));
	vector<vector<double> > p(ny, vector<double>(nx, 0));
	vector<vector<double> > pn(ny, vector<double>(nx, 0));
	vector<vector<double> > b(ny, vector<double>(nx, 0));
	double udiff = 1;	
	int stepcount = 0;
	int num = (ny-2)/(size-1);
	int start = (rank-1)*num+1;
	int end = start+num;
	double buffer[num*2][nx];
	double ucopy[ny][nx];
	double vcopy[ny][nx];
	while(udiff > 0.001){ 	
		for(int i=0;i<un.size();i++){
			for(int j=0;j<un.at(0).size();j++){
				un[i][j] = u[i][j];
				vn[i][j] = v[i][j];
			}
		}

		b = build_up_b(rho,dt,dx,dy,u,v);
		p = pressure_poisson_periodic(p,b,dx,dy,nit);
		int n = u.at(0).size();
		
		//ホストとクライアントに分け、分散処理
		if(rank!=0){
			int l=0;	
			for(int i=start;i<end;i++){
				for(int j=1;j<u.at(0).size()-1;j++){
					buffer[l][j] = (un[i][j]-un[i][j]*dt/dx*(un[i][j]-un[i][j-1])-vn[i][j]*dt/dy*(un[i][j]-un[i-1][j])-dt/(2*rho*dx)*(p[i][j+1]-p[i][j-1])+nu*(dt/pow(dx,2)*(un[i][j+1]-2*un[i][j]+un[i][j-1])+dt/pow(dy,2)*(un[i+1][j]-2*un[i][j]+un[i-1][j]))+F*dt);
				}
				l++;	
			}
	
			n = v.at(0).size();
			for(int i=start;i<end;i++){
				for(int j=1;j<v.at(0).size()-1;j++){
					buffer[l][j] = (vn[i][j]-un[i][j]*dt/dx*(vn[i][j]-vn[i][j-1])-vn[i][j]*dt/dy*(vn[i][j]-vn[i-1][j])-dt/(2*rho*dy)*(p[i+1][j]-p[i-1][j])+nu*(dt/pow(dx,2)*(vn[i][j+1]-2*vn[i][j]+vn[i][j-1])+dt/pow(dy,2)*(vn[i+1][j]-2*vn[i][j]+vn[i-1][j])));
				}
				l++;
			}
    			MPI_Send(buffer, num*2*nx, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    			MPI_Recv(ucopy, ny*nx, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    			MPI_Recv(vcopy, ny*nx, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&udiff,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			for(int i=0;i<ny;i++){
				for(int j=0;j<nx;j++){
					u[i][j] = ucopy[i][j];
					v[i][j] = vcopy[i][j];
				}
			}

		}else{
			for(int i=1;i<u.size()-1;i++){

				u[i][n-1] = (un[i][n-1]-un[i][n-1]*dt/dx*(un[i][n-1]-un[i][n-2])-vn[i][n-1]*dt/dy*(un[i][n-1]-un[i-1][n-1])-dt/(2*rho*dx)*(p[i][0]-p[i][n-2])+nu*(dt/pow(dx,2)*(un[i][0]-2*un[i][n-1]+un[i][n-2])+dt/pow(dy,2)*(un[i+1][n-1]-2*un[i][n-1]+un[i-1][n-1]))+F*dt);
				u[i][0] = (un[i][0]-un[i][0]*dt/dx*(un[i][0]-un[i][n-1])-vn[i][0]*dt/dy*(un[i][0]-un[i-1][0])-dt/(2*rho*dx)*(p[i][1]-p[i][n-1])+nu*(dt/pow(dx,2)*(un[i][1]-2*un[i][0]+un[i][n-1])+dt/pow(dy,2)*(un[i+1][0]-2*un[i][0]+un[i-1][0]))+F*dt);

			}			


			n = v.at(0).size();
			for(int i=1;i<v.size()-1;i++){

				v[i][n-1] = (vn[i][n-1]-un[i][n-1]*dt/dx*(vn[i][n-1]-vn[i][n-2])-vn[i][n-1]*dt/dy*(vn[i][n-1]-vn[i-1][n-1])-dt/(2*rho*dy)*(p[i+1][n-1]-p[i-1][n-1])+nu*(dt/pow(dx,2)*(vn[i][0]-2*vn[i][n-1]+vn[i][n-2])+dt/pow(dy,2)*(vn[i+1][n-1]-2*vn[i][n-1]+vn[i-1][n-1])));
				v[i][0] = (vn[i][0]-un[i][0]*dt/dx*(vn[i][0]-vn[i][n-1])-vn[i][0]*dt/dy*(vn[i][0]-vn[i-1][0])-dt/(2*rho*dy)*(p[i+1][0]-p[i-1][0])+nu*(dt/pow(dx,2)*(vn[i][1]-2*vn[i][0]+vn[i][n-1])+dt/pow(dy,2)*(vn[i+1][0]-2*vn[i][0]+vn[i-1][0])));

			}				
			n = u.at(0).size();
			for(int j=0;j<u.at(0).size();j++){
				u[0][j] = 0;
				u[n-1][j] = 0;
			}
			n = v.at(0).size();
			for(int j=0;j<v.at(0).size();j++){
				v[0][j] = 0;
				v[n-1][j] = 0;
			}
		
	
		//分散処理させた結果をホストで統合
		
			for(int irank=1;irank<size;irank++){
				MPI_Recv(buffer,num*2*nx,MPI_DOUBLE,irank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				for(int i=0;i<num;i++){
					for(int j=1;j<nx-1;j++){
						u[(irank-1)*num+1+i][j] = buffer[i][j];
					}
				}
				for(int i=num;i<num*2;i++){
					for(int j=1;j<nx-1;j++){
						v[(irank-1)*num+1+i-num][j] = buffer[i][j];
					}
				}
			}
			
			double usum = 0;
			double unsum = 0;
			for(int i=0;i<u.size();i++){
				for(int j=0;j<u.at(0).size();j++){
					usum += u[i][j];
					unsum += un[i][j];
				}
			}
		
			udiff = (usum - unsum)/usum;
			stepcount+=1;
			
			for(int i=0;i<ny;i++){
				for(int j=0;j<nx;j++){
					ucopy[i][j] = u[i][j];
					vcopy[i][j] = v[i][j];
				}
			}
			for(int irank=1;irank<size;irank++){
				MPI_Send(ucopy,ny*nx,MPI_DOUBLE,irank,0,MPI_COMM_WORLD);
				MPI_Send(vcopy,ny*nx,MPI_DOUBLE,irank,0,MPI_COMM_WORLD);
				MPI_Send(&udiff,1,MPI_DOUBLE,irank,0,MPI_COMM_WORLD);
			}
		}
	}	
	if(rank==0){cout<<stepcount<<endl;}
	
	string filename = "output.txt";
	ofstream writing_file;
	writing_file.open(filename,ios::out);
	
	for(int i=0;i<ny;i++){
		for(int j=0;j<nx;j++){
			writing_file<<X[i][j]<<endl;
			writing_file<<Y[i][j]<<endl;
			writing_file<<u[i][j]<<endl;
			writing_file<<v[i][j]<<endl;
		}
	}
		
	MPI_Finalize();

}
