#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

#define PI 3.14159265


/*all the unitis here are in Mpc , Solar mass and comovel coordinates*/

/*first lets declare the functions*/

float Deslocation (); /*compute the deslocation of the particles from the zeldovich approximation*/
float VelocityPerturb (); /*compute the comoving perturbation on particles due to the wake*/
void VMaxAndMin (); /*compute VM and Vm from the data*/
void DeslocatedGalaxiesDistribution ();/*Make a dat file with the Galaxies Distribution with a wake*/


/*Now lets define the pointers*/

long long int *countP;               /*count the total number of galaxies in the simulation volume*/
float *zP;     /*redshift*/
float *GmuP;   /*mu is the string tension and G the Newton Constant*/
float *XMP,*XmP,*YMP,*YmP,*ZMP,*ZmP;
float *XRangeWP;   /*is the size of the wake in the X direction*/
float *YRangeWP;   /*is the size of the wake in the Y direction*/
string *filenameP,*filenameOutP;
double *box_sizeP,*n_particlesP;

/*Define the fixed numerical constants*/

float h={0.7};
float OmegaM={0.3};
float OmegaL={0.7};
float clight={9.6E-15};  /*speed of light in Mpc per second units*/
int zi={1000};   /*redshift close to the time of recombination*/
float tzero={4.2E+17}; /*age of the universe today*/
double vSgammaS={clight*(sqrt(3.))/3.}; /*speed of cosmic string times Lorentz factor in Mpc per second units*/
float Hzero={3.2411E-18}; /*hubble constant*/




int main (int argc,char** argv){


srand (time(NULL)); /* initialize random seed: */

/*Define the (initial) values of the pointers*/

float XRangeW;   /*is the size of the simulation box in the X direction*/
float YRangeW;/*is the size of the simulation box in the Y direction*/

float z;     /*redshift*/
float Gmu;   /*mu is the string tension and G the Newton Constant*/
double box_size;
double n_particles;


z=atof(argv[1]);
Gmu=3*pow(10.0,atof(argv[2]));
box_size=atof(argv[7]);
n_particles=atof(argv[8]);

XRangeW=(tzero*clight)/(pow(zi+1,1./2.));
YRangeW=(tzero*vSgammaS)/(pow(zi+1,1./2.));


zP=&z;
GmuP=&Gmu;

XRangeWP=&XRangeW;
YRangeWP=&YRangeW;
box_sizeP=&box_size;
n_particlesP=&n_particles;

printf("The X Size of the wake is %f\t and the Y size is %f\n",*XRangeWP,*YRangeWP );



cout << "z=" << z << endl;
cout << "Wake Thickness" << Deslocation () << endl;
cout << "Velocity perturbation " << VelocityPerturb () << endl;

        string directory_in = argv[3];
        string directory_out = argv[4];
        string filePosName_in = argv[5];
        string filePosName_out = argv[6];

        string filename;
        string filenameOut;

        filename =directory_in + argv[1] + filePosName_in;
        filenameOut = directory_out + argv[1] + filePosName_out;

filenameP=&filename;
filenameOutP=&filenameOut;

cout << (*filenameP).c_str() << endl;

DeslocatedGalaxiesDistribution ();

/*ofstream Info ("/home/acer/Documents/storage/RandomGalaxyDistributionRMR_Info2.dat");

Info << "zi=" << zi << endl;

Info << "z=" << *zP << endl;

Info << "Gmu" << *GmuP << endl;

Info << "Wake Thickness" << Deslocation () << endl;

Info << "XRange= " << *XRangeWP << endl;

Info << "YRange= " << *YRangeWP << endl;

Info.close();
 */

}

float Deslocation () /*compute the deslocation of the galaxies from the zeldovich approximation*/
{
float Desl;

Desl=(float)((12.*3.14)/5.)*(*GmuP)*tzero*vSgammaS*(sqrt(1.+zi))/(1.+(*zP));

return(Desl);
}

float VelocityPerturb () /*compute the comoving velocity perturbation on particles due to the wake*/
{
float VelPert;

VelPert=(float)((8.*3.14)/5.)*(*GmuP)*vSgammaS*(sqrt(1.+zi))*(sqrt(1.+(*zP)));

/* the prefactor (3./2.)*(tzero/sqrt(1+(*zP)))  is introduced to converto to simulation units */

return(VelPert);
}


void DeslocatedGalaxiesDistribution ()/*Make a dat file with the Galaxies Distribution with a wake*/
{

float XM,Xm,YM,Ym,ZM,Zm;
VMaxAndMin ();

XM=*XMP;
Xm=*XmP;
YM=*YMP;
Ym=*YmP;
ZM=*ZMP;
Zm=*ZmP;

cout << "X runs from " << Xm << " to " << XM << " in cell space " << endl;
cout << "Y runs from " << Ym << " to " << YM << " in cell space " << endl;
cout << "Z runs from " << Zm << " to " << ZM << " in cell space " << endl;

int first[12];
float xp[6];

//float LinConv={130./512.}; /*spatial linear convertion factor from simultion to comoving*/
float LinConv=(float){*box_sizeP/(*n_particlesP)}; /*spatial linear convertion factor from simultion to comoving*/
float TimConv={2./(3.*pow(1.+(*zP),2.)*Hzero*pow(OmegaM,1./2.))};/*Convert time in simulation units to seconds*/ 

cout << "X runs from " << Xm*LinConv << " to " << XM*LinConv << " in coordinate space " << endl;
cout << "Y runs from " << Ym*LinConv << " to " << YM*LinConv << " in coordinate space " << endl;
cout << "Z runs from " << Zm*LinConv << " to " << ZM*LinConv << " in coordinate space " << endl;


//ifstream CatalCubep3mInit("/scratch/irulan/camargod/nowake_ran/63.000xv0.dat", ios::binary);

//ofstream CatCubep3mOut("/scratch/irulan/camargod/nowake_ran/63.000xv0_wake.dat");

ifstream CatalCubep3mInit((*filenameP).c_str(), ios::binary);
ofstream CatCubep3mOut ((*filenameOutP).c_str());


/*read and write the header*/

CatalCubep3mInit.read((char*)&first,12*sizeof(int));
CatCubep3mOut.write((char*)&first,12*sizeof(int));

float angle;
angle = -3.14*0./180.;



while (!CatalCubep3mInit.eof())
{

CatalCubep3mInit.read((char*)&xp,6*sizeof(float));

if( CatalCubep3mInit.eof() )
        {
        break;
        }


/*centralize*/



xp[0]=xp[0]-(XM+Xm)/2.;
xp[1]=xp[1]-(YM+Ym)/2.;
xp[2]=xp[2]-(ZM+Zm)/2.;


xp[0] = (xp[0] * (cos ( angle  ))) + (xp[2] * (sin ( angle))) ;
xp[2]=  (xp[2] * (cos ( angle))) - (xp[0] * (sin ( angle  ))) ;

/*move to physical coordinates*/

xp[0]=LinConv*xp[0];
xp[1]=LinConv*xp[1];
xp[2]=LinConv*xp[2];
xp[3]=LinConv*xp[3]/TimConv;
xp[4]=LinConv*xp[4]/TimConv;
xp[5]=LinConv*xp[5]/TimConv;

/*insert the wake*/



	if (xp[0] > (-(*XRangeWP)/2.) && xp[0] < (*XRangeWP)/2.)  {
        if (xp[1] > (-(*YRangeWP)/2.) && xp[1] < (*YRangeWP)/2.)  {
	if (xp[2] > (-(*XRangeWP)/2.) && xp[2] < (*XRangeWP)/2.)  {	
		if (signbit(xp[2]))
     		{
	     	xp[2]=xp[2]+Deslocation ();
	     	xp[5]=xp[5]+VelocityPerturb ();
		}
	  	else
	     	{
	     	xp[2]=xp[2]-Deslocation (); 
	     	xp[5]=xp[5]-VelocityPerturb ();
		}
	}
	}
	}


/*converts back to simulation units*/

xp[0]=xp[0]/LinConv;
xp[1]=xp[1]/LinConv;
xp[2]=xp[2]/LinConv;
xp[3]=TimConv*xp[3]/LinConv;
xp[4]=TimConv*xp[4]/LinConv;
xp[5]=TimConv*xp[5]/LinConv;

/*set initial configuration */

xp[0] = (xp[0] * (cos ( angle  ))) - (xp[2] * (sin ( angle))) ;
xp[2]=  (xp[2] * (cos ( angle))) + (xp[0] * (sin ( angle  ))) ;

xp[0]=xp[0]+(XM+Xm)/2.;
xp[1]=xp[1]+(YM+Ym)/2.;
xp[2]=xp[2]+(ZM+Zm)/2.;


/*write*/

CatCubep3mOut.write((char*)&xp,6*sizeof(float));

}

CatalCubep3mInit.close();
CatCubep3mOut.close();

}


void VMaxAndMin () /*compute VM and Vm from the data*/
{


int first[12];
float xp[6];

//ifstream CatalCubep3mInit("/scratch/irulan/camargod/nowake_ran/63.000xv0.dat", ios::binary);
ifstream CatalCubep3mInit((*filenameP).c_str(), ios::binary);

float PosXM,PosXm,PosYM,PosYm,PosZM,PosZm;

/*read the header*/

CatalCubep3mInit.read((char*)&first,12*sizeof(int));


PosXM=0.;
PosXm=0.;
PosYM=0.;
PosYm=0.;
PosZM=0.;
PosZm=0.;

while (!CatalCubep3mInit.eof()) {

        CatalCubep3mInit.read((char*)&xp,6*sizeof(float));
	if( CatalCubep3mInit.eof() )
        {
        break;
        }


        if(xp[0]>PosXM)
        {
        PosXM=xp[0];
        }
        if(xp[0]<PosXm)
        {
        PosYm=xp[0];
        }

        if(xp[1]>PosYM)
        {
        PosYM=xp[1];
        }
        if(xp[1]<PosYm)
        {
        PosYm=xp[1];
        }

        if(xp[2]>PosZM)
        {
        PosZM=xp[2];
	}
	if(xp[2]<PosZm)
        {
        PosZm=xp[2];
        }
}

CatalCubep3mInit.close();

XMP=&PosXM;
XmP=&PosXm;

YMP=&PosYM;
YmP=&PosYm;

ZMP=&PosZM;
ZmP=&PosZm;
}


