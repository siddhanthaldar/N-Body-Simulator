#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>

using namespace std;

int num_bodies = 1000;

#define n_threads 1
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

	for(i=0; i<num_bodies; i++)
		loc[i] = new float [3];

	i = 0;
	string s, number,tok;
	ifstream f;

	f.open("Trajectory.txt");

	// Ignore first 8 lines
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

// Calculate force on 1 due to 2
float* force(float* loc1, float* loc2) // loc = {x,y,z} 
{
	float* f = new float [3]; //Force
	float dist = 0;
	int i;
	
	for(i=0;i<3;i++)
		dist += (loc2[i]-loc1[i])*(loc2[i]-loc1[i]);

	if(dist==0)
	{
		for(i=0;i<3;i++)
			f[i] = 0;
	}	
	else	
	{
		for(i=0;i<3;i++)
			f[i] = (loc2[i]-loc1[i])/(dist);
	}
	
	return f;
}	
 
void simulate(float** initial_loc)
{
	int i,j,k,l, flag;
	int count = 0;
	float *f, temp;
	ofstream file;
	file.open("output_trajectory.txt");

	//Force
	float** F = new float* [num_bodies];
	for(i=0; i<num_bodies;i++) 
		F[i] = new float [3];

	// Particle velocities
	float** v = new float* [num_bodies];
	for(i=0; i<num_bodies;i++) 
	{
		v[i] = new float [3];
		for(j=0;j<3;j++)
			v[i][j] = 0;
	}	

	// Particle positions -> {num_bodies, time_steps, positions(x,y,z)}
	float*** pos = new float** [num_bodies];
	for(i=0; i<num_bodies;i++)
	{
		pos[i] = new float* [time_steps];
		for(j=0; j<time_steps;j++)
		{
			pos[i][j] = new float [3];
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

	// define cuboid
	int*** cuboid = new int** [length];
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

	// Simulation
	for(i=1;i<time_steps;i++)
	{
		// flag for cuboid to note filled positions
		for(j=0; j<length; j++)
			for(k=0;k<width; k++)
				for(l=0; l<depth; l++)
					cuboid[j][k][l] = -1;

		// Calculate force on each particle
		for(j=0; j<num_bodies; j++)
			for(k=0; k<3; k++)
				F[j][k] = 0;

		for(j=0;j<num_bodies;j++)
			for(k=j+1;k<num_bodies;k++)
			{
				f= force(pos[j][i-1], pos[k][i-1]);
				// Fij = - Fji
				for(l=0;l<3;l++)
				{
					F[j][l] += f[l];
					F[k][l] -= f[l];
				}
			}			

		for(j=0;j<num_bodies;j++)
		{
			flag = 0;
			// Update half step velocity
			for(k=0;k<3;k++)
				v[j][k] += F[j][k]*delta_t/2;  // mass = 1

			// Update particle position
			for(k=0;k<3;k++)
				pos[j][i][k] = pos[j][i-1][k] + v[j][k]*delta_t;
			
			if(pos[j][i][0]<200 && pos[j][i][0]>0 && pos[j][i][1]<100 && pos[j][i][1]>0 && pos[j][i][2]<400 && pos[j][i][2]>0){
				if(cuboid[(int)pos[j][i][1]][(int)pos[j][i][0]][(int)pos[j][i][2]] == -1)  // Position is [width,length,depth]
					cuboid[(int)pos[j][i][1]][(int)pos[j][i][0]][(int)pos[j][i][2]] = j;
				else
					flag = 1;  // Collision occurred
			}

			// Update full step velocity
			for(k=0;k<3;k++)
				v[j][k] += F[j][k]*delta_t/2;  // mass = 1	

			// Detect collision between balls
			if(flag)
			{
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
				if((int)pos[j][i][0] < 0) pos[j][i][0] = 0;
				else if((int)pos[j][i][0]>199) pos[j][i][0] = 199;
			}	
			if((int)pos[j][i][1] <= 0 || (int)pos[j][i][1] >= 99)	// Length
			{
				v[j][1] *= -1.0;
				if((int)pos[j][i][1] < 0) pos[j][i][1] = 0;
				else if((int)pos[j][i][1]>99) pos[j][i][1] = 99;
			}	
			if((int)pos[j][i][2] <= 0 || (int)pos[j][i][2] >= 399)	// Depth
			{
				v[j][2] *= -1.0;
				if((int)pos[j][i][2] < 0) pos[j][i][2] = 0;
				else if((int)pos[j][i][2]>399) pos[j][i][2] = 399;
			}		
		}
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
	float** initial_loc = initial_location();
	cout<<"Started"<<endl;
	simulate(initial_loc);
	cout<<"Num_threads = "<<n_threads<<endl;
	return 0;
}