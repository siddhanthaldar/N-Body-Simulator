#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <omp.h>

using namespace std;

int num_bodies = 1000;

#define n_threads 8
#define length 100
#define width 200
#define depth 400
#define time_steps 720000
float delta_t = 0.01;

// Position is [width,length,depth]

float** initial_location()
{
	float** loc = new float* [num_bodies];
	int i=0,j;

	#pragma omp parallel for shared(loc) private(i) 
	for(i=0; i<num_bodies; i++)
		loc[i] = new float [3];

	i = 0;
	string s, number,tok;
	ifstream f;

	f.open("Trajectory.txt");

	// Ignore first 8 lines
	#pragma omp parallel for private(i) shared(f)
	for(int i=0; i<8; i++)
	{
		f.ignore(1000,'\n');
	}

	i=0;
	while(!f.eof())
	{
		j=0;
		if(i==num_bodies) break;
		getline(f,s);
		stringstream ss(s);
		while(getline(ss,tok, '	')){
			loc[i][j] = atof(tok.c_str());
			j++;
		}
	
		i++;
	}
	f.close();

	return loc;
}

void simulate(float** initial_loc)
{
	int i,j,k,l, flag;
	int count = 0;
	float temp;
	ofstream file;
	file.open("output_trajectory.txt");
	float **F, **v, ***pos;
	int ***cuboid;


	#pragma omp parallel sections private(i,j,k)
	{
        #pragma omp section
		{
			//Force
			F = new float* [num_bodies];
			#pragma omp parallel for shared(F) private(i)
			for(i=0; i<num_bodies;i++) 
				F[i] = new float [3];
		}

        #pragma omp section	
        {
			// Particle velocities
			v = new float* [num_bodies];
			#pragma omp parallel for shared(v) private(i)
			for(i=0; i<num_bodies;i++) 
			{
				v[i] = new float [3];
				for(j=0;j<3;j++)
					v[i][j] = 0;
			}	
        }	

        #pragma omp section
        {
			// Particle positions -> {num_bodies, time_steps, positions(x,y,z)}
			pos = new float** [num_bodies];
			#pragma omp parallel for shared(pos) private(i)
			for(i=0; i<num_bodies;i++)
			{
				pos[i] = new float* [time_steps];
				#pragma omp parallel for shared(pos) private(j)	
				for(j=0; j<time_steps;j++)
				{
					pos[i][j] = new float [3];
				}	
			}	        	
        }

        #pragma omp section
        {
			// define cuboid
			cuboid = new int** [length];
			for(i=0; i<length; i++)
			{
				cuboid[i] = new int* [width];
				for(j=0;j<width; j++)
				{
					cuboid[i][j] = new int [depth];
					for(k=0; k<depth; k++)
						cuboid[i][j][k] = 0;
				}
			}        	
        }
	}

	// Write first time step location to file
	for(i=0; i<num_bodies;i++)
	{
		for(k=0;k<3;k++)
		{
			pos[i][0][k] = initial_loc[i][k];
			file<<pos[i][0][k]<<"	";
		}	
		file<<endl;
	}


	float* f = new float [3]; //Force
	float dist = 0;
	double tot_time = 0;
	double wtime;

	// Simulation
	for(i=1;i<time_steps;i++)
	{
		wtime = omp_get_wtime ( );
		// flag for cuboid to note filled positions
		#pragma omp parallel for collapse(3) shared(cuboid) private(j,k,l)
		for(j=0; j<length; j++)
			for(k=0;k<width; k++)
				for(l=0; l<depth; l++)
					cuboid[j][k][l] = -1;

		// Calculate force on each particle
		#pragma omp parallel for collapse(2) shared(F) private(j,k)
		for(j=0; j<num_bodies; j++)
			for(k=0; k<3; k++)
				F[j][k] = 0;

		// #pragma omp parallel for shared(F,pos) private(j,k,l) lastprivate(f,dist)
		for(j=0;j<num_bodies;j++)
			for(k=j+1;k<num_bodies;k++)
			{
				// Calculate force
				dist=0;
				for(l=0;l<3;l++)
					dist += (pos[k][i-1][l]-pos[j][i-1][l])*(pos[k][i-1][l]-pos[j][i-1][l]);

				if(dist==0)
				{
					for(l=0;l<3;l++)
						f[l] = 0;
				}	
				else	
				{
					for(l=0;l<3;l++)
						f[l] = (pos[k][i-1][l]-pos[j][i-1][l])/(dist);
				}

				// Fij = - Fji
				for(l=0;l<3;l++)
				{
					F[j][l] += f[l];
					F[k][l] -= f[l];
				}
			}			

		#pragma omp parallel for shared(F,v,pos,cuboid, delta_t) private(j,k,l,flag,temp)
		for(j=0;j<num_bodies;j++)
		{
			flag = 0;
			// Update half step velocity
			#pragma omp parallel for shared(v,j,delta_t) private(k)
			for(k=0;k<3;k++)
				v[j][k] += F[j][k]*delta_t/2;  // mass = 1

			// Update particle position
			// #pragma omp parallel for shared(pos,j,i,delta_t) private(k)
			for(k=0;k<3;k++)
				pos[j][i][k] = pos[j][i-1][k] + v[j][k]*delta_t;
			
			if(pos[j][i][0]<200 && pos[j][i][0]>0 && pos[j][i][1]<100 && pos[j][i][1]>0 && pos[j][i][2]<400 && pos[j][i][2]>0){
				if(cuboid[(int)pos[j][i][1]][(int)pos[j][i][0]][(int)pos[j][i][2]] == -1)  // Position is [width,length,depth]
					cuboid[(int)pos[j][i][1]][(int)pos[j][i][0]][(int)pos[j][i][2]] = j;
				else
					flag = 1;  // Collision occurred
			}

			// Update full step velocity
			#pragma omp parallel for shared(v,j,delta_t) private(k)
			for(k=0;k<3;k++)
				v[j][k] += F[j][k]*delta_t/2;  // mass = 1	

			// Detect collision between balls
			if(flag)
			{
				#pragma omp parallel for shared(v,cuboid,pos,j,i,delta_t) private(temp,k)
				for(k=0;k<3;k++)
				{
					temp = v[j][k];
					v[j][k] = v[cuboid[(int)pos[j][i][1]][(int)pos[j][i][0]][(int)pos[j][i][2]]][k];
					v[cuboid[(int)pos[j][i][1]][(int)pos[j][i][0]][(int)pos[j][i][2]]][k]= temp;  
				}	
				continue;
			}

			// Detect Collision with wall 
					if((int)pos[j][i][0] <= 0 || (int)pos[j][i][0] >= 199)	// Width
					{
						v[j][0] *= -1.0;
						if((int)pos[j][i][0] < 0) pos[j][i][0] *= -1;
						else if((int)pos[j][i][0]>199) pos[j][i][0] -= 199;
					}	

					if((int)pos[j][i][1] <= 0 || (int)pos[j][i][1] >= 99)	// Length
					{
						v[j][1] *= -1.0;
						if((int)pos[j][i][1] < 0) pos[j][i][1] *= -1;
						else if((int)pos[j][i][1]>99) pos[j][i][1] -= 99;
					}						

					if((int)pos[j][i][2] <= 0 || (int)pos[j][i][2] >= 399)	// Depth
					{
						v[j][2] *= -1.0;
						if((int)pos[j][i][2] < 0) pos[j][i][2] *= -1;
						else if((int)pos[j][i][2]>399) pos[j][i][2] -= 399;
					}		
	
		}

		wtime = omp_get_wtime ( ) - wtime;
		tot_time += wtime;
		cout<<"Step : "<<i<<"	Time = "<<wtime<<"	Total time = "<<tot_time<<endl;

		// write pos to file
		if(i%100==0)
		{
			cout<<"Time step : "<<i<<endl;
			for(j=0;j<num_bodies;j++)
			{
				for(k=0;k<3;k++)
					file<<pos[j][i][k]<<"	";
				file<<endl;
			}
		}		
	}

	file.close();
}

int main()
{
	omp_set_num_threads(n_threads);
	// double wtime;
	// wtime = omp_get_wtime ( );
	float** initial_loc = initial_location();
	simulate(initial_loc);
	// wtime = omp_get_wtime ( ) - wtime;
	// cout<<"Num_threads = "<<n_threads<<"	Time = "<<wtime<<endl;
	return 0;
}