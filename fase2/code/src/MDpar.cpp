/*
 MD.c - a simple molecular dynamics program for simulating real gas properties of Lennard-Jones particles.
 
 Copyright (C) 2016  Jonathan J. Foley IV, Chelsea Sweet, Oyewumi Akinfenwa
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
 Electronic Contact:  foleyj10@wpunj.edu
 Mail Contact:   Prof. Jonathan Foley
 Department of Chemistry, William Paterson University
 300 Pompton Road
 Wayne NJ 07470
 
 */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<vector>
#include<omp.h>


// Number of particles
const int N = 10*500;

//  Lennard-Jones parameters in natural units!
double sigma = 1.;
double epsilon = 1.;
double m = 1.;
double kB = 1.;

double NA = 6.022140857e23;
double kBSI = 1.38064852e-23;  // m^2*kg/(s^2*K)

//  Size of box, which will be specified in natural units
double L;

// Global variables to store kinetic and potential values, used in MSV_Kinetic and VV_Pot
double Global_KE;
double Global_Pot;

//  Initial Temperature in Natural Units
double Tinit;
const int MAXPART=5001;
const int MAXPART3 = MAXPART * 3;

// Changed matrix into a vector
//  Position
double r[MAXPART3];
//  Velocity
double v[MAXPART3];
//  Acceleration
double a[MAXPART3];
//  Force
double F[MAXPART3];

// atom type
char atype[10];
//  Function prototypes
//  initialize positions on simple cubic lattice, also calls function to initialize velocities
void initialize();  
//  update positions and velocities using Velocity Verlet algorithm 
//  print particle coordinates to file for rendering via VMD or other animation software
//  return 'instantaneous pressure'
double VelocityVerlet(double dt, int iter, FILE *fp);  
//  Compute Force using F = -dV/dr
//  solve F = ma for use in Velocity Verlet
void computeAccelerations();
//  Numerical Recipes function for generation gaussian distribution
double gaussdist();
//  Initialize velocities according to user-supplied initial Temperature (Tinit)
void initializeVelocities();
//  Compute total potential energy from particle coordinates
double Potential();
double Potential2();
//  Compute mean squared velocity from particle velocities
double MeanSquaredVelocity();
//  Compute total kinetic energy from particle mass and velocities
double Kinetic();

// Junction of ComputeAccelerations() and Potential()
double computeAccelerationsAndPot();

// Junction of MeanSquaredVelocity() and Kinetic()
double MSV_Kinetic();

// Junction of VelocityVerlet() and Potential()
double VV_Pot(double, int, FILE*);

// self explanatory
const char * OUTPUT_FOLDER = "output/";

int main()
{
    
    //  variable delcarations
    int i;
    double dt, Vol, Temp, Press, Pavg, Tavg, rho;
    double VolFac, TempFac, PressFac, timefac;
    double KE, PE, mvs, gc, Z;
    char prefix[1000], tfn[1000], ofn[1000], afn[1000];
    FILE *tfp, *ofp, *afp;
    
    
    printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("                  WELCOME TO WILLY P CHEM MD!\n");
    printf("  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("\n  ENTER A TITLE FOR YOUR CALCULATION!\n");
    
    if (scanf("%s", prefix) != 1) {
        fprintf(stderr, "Error reading prefix.\n");
        exit(1);
    }

    strcpy(tfn,OUTPUT_FOLDER);
    strcat(tfn,prefix);
    strcat(tfn,"_traj.xyz");

    strcpy(ofn,OUTPUT_FOLDER);
    strcat(ofn,prefix);
    strcat(ofn,"_output.txt");

    strcpy(afn,OUTPUT_FOLDER);
    strcat(afn,prefix);
    strcat(afn,"_average.txt");
    
    printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("                  TITLE ENTERED AS '%s'\n",prefix);
    printf("  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    
    /*     Table of values for Argon relating natural units to SI units:
     *     These are derived from Lennard-Jones parameters from the article
     *     "Liquid argon: Monte carlo and molecular dynamics calculations"
     *     J.A. Barker , R.A. Fisher & R.O. Watts
     *     Mol. Phys., Vol. 21, 657-673 (1971)
     *
     *     mass:     6.633e-26 kg          = one natural unit of mass for argon, by definition
     *     energy:   1.96183e-21 J      = one natural unit of energy for argon, directly from L-J parameters
     *     length:   3.3605e-10  m         = one natural unit of length for argon, directly from L-J parameters
     *     volume:   3.79499-29 m^3        = one natural unit of volume for argon, by length^3
     *     time:     1.951e-12 s           = one natural unit of time for argon, by length*sqrt(mass/energy)
     ***************************************************************************************/
    
    //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //  Edit these factors to be computed in terms of basic properties in natural units of
    //  the gas being simulated
    
    
    printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("  WHICH NOBLE GAS WOULD YOU LIKE TO SIMULATE? (DEFAULT IS ARGON)\n");
    printf("\n  FOR HELIUM,  TYPE 'He' THEN PRESS 'return' TO CONTINUE\n");
    printf("  FOR NEON,    TYPE 'Ne' THEN PRESS 'return' TO CONTINUE\n");
    printf("  FOR ARGON,   TYPE 'Ar' THEN PRESS 'return' TO CONTINUE\n");
    printf("  FOR KRYPTON, TYPE 'Kr' THEN PRESS 'return' TO CONTINUE\n");
    printf("  FOR XENON,   TYPE 'Xe' THEN PRESS 'return' TO CONTINUE\n");
    printf("  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    if (scanf("%s", atype) != 1) {
        fprintf(stderr, "Error reading atype.\n");
        exit(1);
    }

    
    if (strcmp(atype,"He")==0) {
        VolFac = 1.8399744000000005e-29;
        PressFac = 8152287.336171632;
        TempFac = 10.864459551225972;
        timefac = 1.7572698825166272e-12;
    }
    else if (strcmp(atype,"Ne")==0) {
        VolFac = 2.0570823999999997e-29;
        PressFac = 27223022.27659913;
        TempFac = 40.560648991243625;
        timefac = 2.1192341945685407e-12;
    }
    else if (strcmp(atype,"Ar")==0) {
        VolFac = 3.7949992920124995e-29;
        PressFac = 51695201.06691862;
        TempFac = 142.0950000000000;
        timefac = 2.09618e-12;
        //strcpy(atype,"Ar");
    }
    else if (strcmp(atype,"Kr")==0) {
        VolFac = 4.5882712000000004e-29;
        PressFac = 59935428.40275003;
        TempFac = 199.1817584391428;
        timefac = 8.051563913585078e-13;
    }
    else if (strcmp(atype,"Xe")==0) {
        VolFac = 5.4872e-29;
        PressFac = 70527773.72794868;
        TempFac = 280.30305642163006;
        timefac = 9.018957925790732e-13;
    }
    else {
        VolFac = 3.7949992920124995e-29;
        PressFac = 51695201.06691862;
        TempFac = 142.0950000000000;
        timefac = 2.09618e-12;
        strcpy(atype,"Ar");
    }
    printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("\n                     YOU ARE SIMULATING %s GAS! \n",atype);
    printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    
    printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("\n  YOU WILL NOW ENTER A FEW SIMULATION PARAMETERS\n");
    printf("  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("\n\n  ENTER THE INTIAL TEMPERATURE OF YOUR GAS IN KELVIN\n");
    if (scanf("%lf", &Tinit) != 1) {
        fprintf(stderr, "Error reading Tinit.\n");
        exit(1);
    }

    // Make sure temperature is a positive number!
    if (Tinit<0.) {
        printf("\n  !!!!! ABSOLUTE TEMPERATURE MUST BE A POSITIVE NUMBER!  PLEASE TRY AGAIN WITH A POSITIVE TEMPERATURE!!!\n");
        exit(0);
    }
    // Convert initial temperature from kelvin to natural units
    Tinit /= TempFac;
    
    
    printf("\n\n  ENTER THE NUMBER DENSITY IN moles/m^3\n");
    printf("  FOR REFERENCE, NUMBER DENSITY OF AN IDEAL GAS AT STP IS ABOUT 40 moles/m^3\n");
    printf("  NUMBER DENSITY OF LIQUID ARGON AT 1 ATM AND 87 K IS ABOUT 35000 moles/m^3\n");
    
    if (scanf("%lf", &rho) != 1) {
        fprintf(stderr, "Error reading rho.\n");
        exit(1);
    }

    Vol = N/(rho*NA);
    
    Vol /= VolFac;
    
    //  Limiting N to MAXPART for practical reasons
    if (N>=MAXPART) {
        
        printf("\n\n\n  MAXIMUM NUMBER OF PARTICLES IS %i\n\n  PLEASE ADJUST YOUR INPUT FILE ACCORDINGLY \n\n", MAXPART);
        exit(0);
        
    }
    //  Check to see if the volume makes sense - is it too small?
    //  Remember VDW radius of the particles is 1 natural unit of length
    //  and volume = L*L*L, so if V = N*L*L*L = N, then all the particles
    //  will be initialized with an interparticle separation equal to 2xVDW radius
    if (Vol<N) {
        
        printf("\n\n\n  YOUR DENSITY IS VERY HIGH!\n\n");
        printf("  THE NUMBER OF PARTICLES IS %i AND THE AVAILABLE VOLUME IS %f NATURAL UNITS\n",N,Vol);
        printf("  SIMULATIONS WITH DENSITY GREATER THAN 1 PARTCICLE/(1 Natural Unit of Volume) MAY DIVERGE\n");
        printf("  PLEASE ADJUST YOUR INPUT FILE ACCORDINGLY AND RETRY\n\n");
        exit(0);
    }
    // Vol = L*L*L;
    // Length of the box in natural units:
    L = pow(Vol,(1./3));
    
    //  Files that we can write different quantities to
    tfp = fopen(tfn,"w");     //  The MD trajectory, coordinates of every particle at each timestep
    ofp = fopen(ofn,"w");     //  Output of other quantities (T, P, gc, etc) at every timestep
    afp = fopen(afn,"w");    //  Average T, P, gc, etc from the simulation
    
    int NumTime;
    if (strcmp(atype,"He")==0) {
        
        // dt in natural units of time s.t. in SI it is 5 f.s. for all other gasses
        dt = 0.2e-14/timefac;
        //  We will run the simulation for NumTime timesteps.
        //  The total time will be NumTime*dt in natural units
        //  And NumTime*dt multiplied by the appropriate conversion factor for time in seconds
        NumTime=50000;
    }
    else {
        dt = 0.5e-14/timefac;
        NumTime=200;
        
    }
    
    //  Put all the atoms in simple crystal lattice and give them random velocities
    //  that corresponds to the initial temperature we have specified
    initialize();
    
    //  Based on their positions, calculate the ininial intermolecular forces
    //  The accellerations of each particle will be defined from the forces and their
    //  mass, and this will allow us to update their positions via Newton's law
    computeAccelerations();
    
    
    // Print number of particles to the trajectory file
    fprintf(tfp,"%i\n",N);
    
    //  We want to calculate the average Temperature and Pressure for the simulation
    //  The variables need to be set to zero initially
    Pavg = 0;
    Tavg = 0;
    
    
    int tenp = floor(NumTime/10);
    fprintf(ofp,"  time (s)              T(t) (K)              P(t) (Pa)           Kinetic En. (n.u.)     Potential En. (n.u.) Total En. (n.u.)\n");
    printf("  PERCENTAGE OF CALCULATION COMPLETE:\n  [");
    for (i=0; i<NumTime+1; i++) {
        
        //  This just prints updates on progress of the calculation for the users convenience
        if (i==tenp) printf(" 10 |");
        else if (i==2*tenp) printf(" 20 |");
        else if (i==3*tenp) printf(" 30 |");
        else if (i==4*tenp) printf(" 40 |");
        else if (i==5*tenp) printf(" 50 |");
        else if (i==6*tenp) printf(" 60 |");
        else if (i==7*tenp) printf(" 70 |");
        else if (i==8*tenp) printf(" 80 |");
        else if (i==9*tenp) printf(" 90 |");
        else if (i==10*tenp) printf(" 100 ]\n");
        fflush(stdout);
        
        
        // This updates the positions and velocities using Newton's Laws
        // Also computes the Pressure as the sum of momentum changes from wall collisions / timestep
        // which is a Kinetic Theory of gasses concept of Pressure
        
        /*
        // VERSÃO SEPARADA
        Press = VelocityVerlet(dt, i+1, tfp);
        Press *= PressFac;
        PE = Potential();
        //PE = Potential2();
        */
        
        
        // VERSÃO JUNTA
        Press = VV_Pot(dt, i+1, tfp);
        Press *= PressFac;
        // Updated the potential value
        PE = Global_Pot;
        
 
        //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //  Now we would like to calculate somethings about the system:
        //  Instantaneous mean velocity squared, Temperature, Pressure
        //  Potential, and Kinetic Energy
        //  We would also like to use the IGL to try to see if we can extract the gas constant

        mvs = MSV_Kinetic();
        // Updated the kinetic value
        KE = Global_KE;

        // Temperature from Kinetic Theory
        Temp = m*mvs/(3*kB) * TempFac;
        
        // Instantaneous gas constant and compressibility - not well defined because
        // pressure may be zero in some instances because there will be zero wall collisions,
        // pressure may be very high in some instances because there will be a number of collisions
        gc = NA*Press*(Vol*VolFac)/(N*Temp);
        Z  = Press*(Vol*VolFac)/(N*kBSI*Temp);
        
        Tavg += Temp;
        Pavg += Press;
        
        fprintf(ofp,"  %8.4e  %20.12g  %20.12g %20.12g  %20.12g  %20.12g \n",i*dt*timefac,Temp,Press,KE, PE, KE+PE);
    }
    
    // Because we have calculated the instantaneous temperature and pressure,
    // we can take the average over the whole simulation here
    Pavg /= NumTime;
    Tavg /= NumTime;
    Z = Pavg*(Vol*VolFac)/(N*kBSI*Tavg);
    gc = NA*Pavg*(Vol*VolFac)/(N*Tavg);
    fprintf(afp,"     Total Time (s)                 T (K)                         P (Pa)                PV/nT (J/(mol K))                  Z                        V (m^3)                       N\n");
    fprintf(afp," -----------------------   ----------------------       -------------------------   -------------------------   ------------------------   -----------------------   ---------------------------\n");
    fprintf(afp,"  %8.12e       %15.12g              %15.12g       %10.12g              %10.12g             %10.12e          %i\n",i*dt*timefac,Tavg,Pavg,gc,Z,Vol*VolFac,N);
    
    printf("\n  TO ANIMATE YOUR SIMULATION, OPEN THE FILE \n  '%s' WITH VMD AFTER THE SIMULATION COMPLETES\n",tfn);
    printf("\n  TO ANALYZE INSTANTANEOUS DATA ABOUT YOUR MOLECULE, OPEN THE FILE \n  '%s' WITH YOUR FAVORITE TEXT EDITOR OR IMPORT THE DATA INTO EXCEL\n",ofn);
    printf("\n  THE FOLLOWING THERMODYNAMIC AVERAGES WILL BE COMPUTED AND WRITTEN TO THE FILE  \n  '%s':\n",afn);
    printf("\n  AVERAGE TEMPERATURE (K):                 %15.12f\n",Tavg);
    printf("\n  AVERAGE PRESSURE  (Pa):                  %15.12f\n",Pavg);
    printf("\n  PV/nT (J * mol^-1 K^-1):                 %15.12f\n",gc);
    printf("\n  PERCENT ERROR of pV/nT AND GAS CONSTANT: %15.12f\n",100*fabs(gc-8.3144598)/8.3144598); // NOTA: molar gas constant,
    printf("\n  THE COMPRESSIBILITY (unitless):          %15.12f \n",Z);
    printf("\n  TOTAL VOLUME (m^3):                      %10.12e \n",Vol*VolFac);
    printf("\n  NUMBER OF PARTICLES (unitless):          %i \n", N);
    
    fclose(tfp);
    fclose(ofp);
    fclose(afp);
    
    return 0;
}


void initialize() {
    int n, p, i, j, k;
    double pos;
    
    // Number of atoms in each direction
    n = int(ceil(pow(N, 1.0/3)));
    
    //  spacing between atoms along a given direction
    pos = L / n;
    
    //  index for number of particles assigned positions
    p = 0;
    //  initialize positions
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            for (k=0; k<n; k++) {
                if (p<N*3) {
                    r[p] = (i + 0.5)*pos;
                    r[p+1] = (j + 0.5)*pos;
                    r[p+2] = (k + 0.5)*pos;
                }
                p+=3;
            }
        }
    }
    
    // Call function to initialize velocities
    initializeVelocities();
    
    /***********************************************
    // Uncomment if you want to see what the initial positions and velocities are
    printf("  Printing initial positions!\n");
    for (i=0; i<N; i++) {
    printf("  %6.3e  %6.3e  %6.3e\n",r[i][0],r[i][1],r[i][2]);
    }
    
    printf("  Printing initial velocities!\n");
    for (i=0; i<N; i++) {
    printf("  %6.3e  %6.3e  %6.3e\n",v[i][0],v[i][1],v[i][2]);
    }
    */
}   

// Junction of MeanSquaredVelocity() and Kinetic()
double MSV_Kinetic(){
    
    double vx2 = 0;
    double vy2 = 0;
    double vz2 = 0;
    double v22;

    double v0, v1, v2;

    double msv, kin;
    
    kin =0.;
    for (int i=0; i<N*3; i+=3) {
        v0 = v[i]*v[i];
        v1 = v[i+1]*v[i+1];
        v2 = v[i+2]*v[i+2];
        
        vx2 += v0;
        vy2 += v1;
        vz2 += v2;

        v22 = v0 + v1 + v2;
        
        kin += m*v22/2.;
    }
    msv = (vx2+vy2+vz2)/N;
    
    //printf("  Average of x-component of velocity squared is %f\n",v2);
    // Updated the kinetic value
    Global_KE = kin;
    return msv;
}

// paralelizada, sem o if
double Potential() {
    double Pot=0;
    int i, j;
    Pot=0.;
    
    #pragma omp parallel for reduction(+:Pot) schedule(dynamic)
    for (i=0; i<N*3; i+=3) {
        for (j = i + 3; j < N*3; j+=3) {
            // Calculate differences in positions
            double dx = r[i] - r[j];
            double dy = r[i+1] - r[j+1];
            double dz = r[i+2] - r[j+2];

            // Calculate the squared distance
            double rSqd = dx * dx + dy * dy + dz * dz;

            // Calculate inverse powers of rSqd
            double inv_rSqd = 1.0 / rSqd;
            double inv_rSqd_2 = inv_rSqd * inv_rSqd;
            double inv_rSqd_6 = inv_rSqd_2 * inv_rSqd_2 * inv_rSqd_2;
            double inv_rSqd_3 = inv_rSqd * inv_rSqd * inv_rSqd;

            Pot += (inv_rSqd_6 - inv_rSqd_3)*2;
            

        }
    }
    return 4.0 * epsilon * Pot;
}

// PARALELIZA CORRETAMENTE !!!
void computeAccelerations() {
    int i, j;
    double cache[3] = {0};

    // Loop over all distinct pairs i,j
    #pragma omp parallel for private(cache) reduction(+:a) schedule(dynamic)
    for (i = 0; i < N * 3; i += 3) {
        cache[0] = 0;
        cache[1] = 0;
        cache[2] = 0;
        for (j = i + 3; j < N * 3; j += 3) {
            // Calculate differences in positions
            double dx = r[i] - r[j];
            double dy = r[i+1] - r[j+1];
            double dz = r[i+2] - r[j+2];

            // Calculate the squared distance
            double rSqd = dx * dx + dy * dy + dz * dz;

            // Calculate inverse powers of rSqd
            double inv_rSqd = 1.0 / rSqd;
            double inv_rSqd_2 = inv_rSqd * inv_rSqd;
            double inv_rSqd_4 = inv_rSqd_2 * inv_rSqd_2;
            double inv_rSqd_7 = inv_rSqd * inv_rSqd_2 * inv_rSqd_4;

            // Calculate the force magnitude
            double f = (48 * inv_rSqd_7 - 24 * inv_rSqd_4);

            cache[0] += f * dx;
            cache[1] += f * dy;
            cache[2] += f * dz;

            a[j] -= f*dx;
            a[j+1] -= f*dy;
            a[j+2] -= f*dz;
        }
        // Update accelerations
        a[i] += cache[0];
        a[i+1] += cache[1];
        a[i+2] += cache[2];
    }
}

double computeAccelerationsAndPot(){
    // Potential
    double Pot;
    double cache[3] = {0};
    int i, j;
    Pot = 0.0;
    
    // loop over all distinct pairs i,j
    #pragma omp parallel for private(cache) reduction(+:Pot) reduction(+:a) schedule(dynamic)
    for (i = 0; i < (N-1) * 3; i += 3) {
        cache[0] = 0;
        cache[1] = 0;
        cache[2] = 0;
        for (j = i + 3; j < N * 3; j += 3) {
            // Calculate differences in positions
            double dx = r[i] - r[j];
            double dy = r[i+1] - r[j+1];
            double dz = r[i+2] - r[j+2];

            // Calculate the squared distance
            double rSqd = dx * dx + dy * dy + dz * dz;

            // Calculate inverse powers of rSqd
            double inv_rSqd = 1.0 / rSqd;
            double inv_rSqd_2 = inv_rSqd * inv_rSqd;
            double inv_rSqd_3 = inv_rSqd * inv_rSqd * inv_rSqd;
            double inv_rSqd_4 = inv_rSqd_2 * inv_rSqd_2;
            double inv_rSqd_6 = inv_rSqd_2 * inv_rSqd_2 * inv_rSqd_2;
            double inv_rSqd_7 = inv_rSqd * inv_rSqd_2 * inv_rSqd_4;

            // Calculate the force magnitude
            double f = (48 * inv_rSqd_7 - 24 * inv_rSqd_4);

            cache[0]+=f*dx;
            cache[1]+=f*dy;
            cache[2]+=f*dz;

            a[j] -= f * dx;
            a[j+1] -= f * dy;
            a[j+2] -= f * dz;

            // ---------------------- begin pot
            Pot += (inv_rSqd_6 - inv_rSqd_3) * 2;
            // ---------------------- end pot
        }
        // update accelerations of i
        a[i] += cache[0];
        a[i+1] += cache[1];
        a[i+2] += cache[2];
    }
    return Pot;
}

double VV_Pot(double dt, int iter, FILE *fp){
    // Velocity Verlet
    double psum = 0.;
    double dt2 = dt*dt;
    double dt5 = dt * 0.5;
    double aij5[3];
    double m2 = 2*m;
    int i;

    double Pot;

    //#pragma omp parallel for private(aij5)
    for (i=0; i<N*3; i+=3) {
        aij5[0] =  0.5*a[i];
        aij5[1] =  0.5*a[i+1];
        aij5[2] =  0.5*a[i+2];

        r[i] += v[i]*dt + aij5[0]*dt2;
        r[i+1] += v[i+1]*dt + aij5[1]*dt2;
        r[i+2] += v[i+2]*dt + aij5[2]*dt2;

        v[i] += aij5[0]*dt;
        v[i+1] += aij5[1]*dt;
        v[i+2] += aij5[2]*dt;

        a[i] = 0;
        a[i+1] = 0;
        a[i+2]= 0;
    }

    //  Update accellerations from updated positions
    Pot = computeAccelerationsAndPot();

    //  Update velocity with updated acceleration
    //#pragma omp parallel for reduction(+:psum)
    for (i=0; i<N*3; i+=3) {
        v[i] += a[i]*dt5;
        v[i+1] += a[i+1]*dt5;
        v[i+2] += a[i+2]*dt5;

        // Elastic walls
        if (r[i] < 0. || r[i] >= L) {
            v[i] *= -1.;
            psum += m2*fabs(v[i])/dt;// contribution to pressure from "left" and "right" walls
        }
        if (r[i+1] < 0. || r[i+1] >= L) {
            v[i+1] *= -1.;
            psum += m2*fabs(v[i+1])/dt;// contribution to pressure from "left" and "right" walls
        }
        if (r[i+2] < 0. || r[i+2] >= L) {
            v[i+2] *= -1.;
            psum += m2*fabs(v[i+2])/dt;// contribution to pressure from "left" and "right" walls
        }
    }
    
    double vv = psum/(6*L*L);;
    Global_Pot = 4.0 * epsilon * Pot;
    return vv;
}


// -----------------------------------------------------------------------------------------
// returns sum of dv/dt*m/A (aka Pressure) from elastic collisions with walls
double VelocityVerlet(double dt, int iter, FILE *fp) {
    int i;
    
    double psum = 0.;

    double dt2 = dt*dt;
    double dt5 = dt * 0.5;

    double aij5[3];

    double m2 = 2*m;
    
    #pragma omp parallel for private(aij5) // TODO
    for (i=0; i<N*3; i+=3) {
        aij5[0] =  0.5*a[i];
        aij5[1] =  0.5*a[i+1];
        aij5[2] =  0.5*a[i+2];

        r[i] += v[i]*dt + aij5[0]*dt2;
        r[i+1] += v[i+1]*dt + aij5[1]*dt2;
        r[i+2] += v[i+2]*dt + aij5[2]*dt2;

        v[i] += aij5[0]*dt;
        v[i+1] += aij5[1]*dt;
        v[i+2] += aij5[2]*dt;

        a[i] = 0;
        a[i+1] = 0;
        a[i+2]= 0;
    }

    //  Update accellerations from updated positions
    computeAccelerations();

    //  Update velocity with updated acceleration
    #pragma omp parallel for reduction(+:psum) // TODO
    for (i=0; i<N*3; i+=3) {
        v[i] += a[i]*dt5;
        v[i+1] += a[i+1]*dt5;
        v[i+2] += a[i+2]*dt5;

        // Elastic walls
        if (r[i] < 0. || r[i] >= L) {
            v[i] *= -1.;
            psum += m2*fabs(v[i])/dt;// contribution to pressure from "left" and "right" walls
        }
        if (r[i+1] < 0. || r[i+1] >= L) {
            v[i+1] *= -1.;
            psum += m2*fabs(v[i+1])/dt;// contribution to pressure from "left" and "right" walls
        }
        if (r[i+2] < 0. || r[i+2] >= L) {
            v[i+2] *= -1.;
            psum += m2*fabs(v[i+2])/dt;// contribution to pressure from "left" and "right" walls
        }
    }
    return psum/(6*L*L);
}

void initializeVelocities() {
    
    int i;
    
    for (i=0; i<N*3; i+=3) {
        //  Pull a number from a Gaussian Distribution
        v[i] = gaussdist();
        v[i+1] = gaussdist();
        v[i+2] = gaussdist();
    }
    
    double vCM[3] = {0, 0, 0};
    
    for (i=0; i<N*3; i+=3) {
        vCM[0] += m*v[i];
        vCM[1] += m*v[i+1];
        vCM[2] += m*v[i+2];
    }
    
    for (i=0; i<3; i++) vCM[i] /= N*m;
    
    for (i=0; i<N*3; i+=3) {
        v[i] -= vCM[0];
        v[i+1] -= vCM[1];
        v[i+2] -= vCM[2];
    }
    
    double vSqdSum, lambda;
    vSqdSum=0.;
    for (i=0; i<N*3; i+=3) {
        vSqdSum += v[i]*v[i];
        vSqdSum += v[i+1]*v[i+1];
        vSqdSum += v[i+2]*v[i+2];
    }
    
    lambda = sqrt( 3*(N-1)*Tinit/vSqdSum);
    
    for (i=0; i<N*3; i+=3) {
        v[i] *= lambda;
        v[i+1] *= lambda;
        v[i+2] *= lambda;
    }
}


//  Numerical recipes Gaussian distribution number generator
double gaussdist() {
    static bool available = false;
    static double gset;
    double fac, rsq, v1, v2;
    if (!available) {
        do {
            v1 = 2.0 * rand() / double(RAND_MAX) - 1.0;
            v2 = 2.0 * rand() / double(RAND_MAX) - 1.0;
            rsq = v1 * v1 + v2 * v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        
        fac = sqrt(-2.0 * log(rsq) / rsq);
        gset = v1 * fac;
        available = true;
        
        return v2*fac;
    } else {
        
        available = false;
        return gset;
    }
}
