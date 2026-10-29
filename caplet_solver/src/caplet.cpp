/*
Created : 2010-07-28
Author  : Yu-Chung Hsiao
Email   : project.caplet@gmail.com
*/

/*
This file is part of CAPLET.

CAPLET is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

CAPLET is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CAPLET.  If not, see <http://www.gnu.org/licenses/>.
*/



#include "caplet_const.h"
#include "caplet_parameter.h"
#include "caplet.h"
#include "caplet_blas.h"
#include "caplet_widgets.h"
#include "caplet_elem.h"
#include "caplet_int.h"
#include "caplet_gauss.h"

#include "mpi.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <list>
#include <cmath>
#include <cstdlib>
#include <iomanip>


namespace caplet{

using namespace std;

//**** Private constants
float Caplet::epsilon0 	= 8.8541878176e-12f;
float Caplet::pi		= 3.1415926535f;

//__________________________________________________________
//*
//* Ctor and Dtor
//*

Caplet::Caplet()
    : isLoaded(false), isSolved(false),
      nDiels(0), DielPanels(0), nDielPanels(0), nDielCoefs(0),
      flagMergeProjection1_0(true){
}


Caplet::~Caplet(){
    this->clear();
}


void Caplet::clear(){
    if ( this->isLoaded == true){
	    this->isLoaded = false;
        this->isInstantiable = false;

        delete[] this->nWirePanels;
        delete[] this->nWireCoefs;
        if (this->nDiels > 0) {
            delete[] this->nDielPanels;
            delete[] this->nDielCoefs;
        }
        this->isDielectric.clear();
        this->materialA.clear();
        this->materialB.clear();
        this->nDiels = 0;
        this->DielPanels = 0;

        delete[] this->panels;

        delete[] this->dirs;
        delete[] this->areas;
        delete[] this->indexIncrements;
        delete[] this->basisTypes;
        delete[] this->basisDirs;
        delete[] this->basisZs;
        delete[] this->basisShifts;

		if ( this->isSolved == true ){
			this->isSolved = false;

            //* single precision version
            switch(this->mode){
            case FAST_GALERKIN:
                delete[] this->P;
                delete[] this->rhs;
                delete[] this->coefs;

                delete[] this->Cmat;
                break;
            case DOUBLE_GALERKIN:
            case DOUBLE_COLLOCATION:
                delete[] this->dP;
                delete[] this->drhs;
                delete[] this->dcoefs;

                delete[] this->dCmat;
                break;
            default:
                ;
            }
        }
	}
    this->isInstantiable = false;
}


//__________________________________________________________
//*
//* Versioned algorithms
//*


//__________________________________________________________
//*
//* Public functions
//*

void Caplet::loadFastcapFile(const std::string filename){

    //* Clean up
    this->clear();

    ifstream ifile(filename.c_str());
    //* Check if open successfully
	if (!ifile){
        cerr << "ERROR: Cannot open the file: " << filename << endl;
		return;
	}
    
    //* Initialize
	this->nPanels = 0;
	this->nWires = 0;

    string		stringTemp;
    string		lineTemp;
    string      currentConductorName = "";
    char        charTemp;
    list<int>	nWirePanelsList;
    int			nWirePanelTemp = 0;

    string condString = "";
    string dielString = "";

    /*Dielectrics*/
    string currentDielectricName = "";
    int nDielPanelTemp = 0;
    std::vector<int> nDielPanelsList;
    this->nDiels = 0;
    double epsA, epsB;

    this->nTotalDielPanels=0;
    this->nTotalWirePanels=0;




	float xCoord[4];
	float yCoord[4];
	float zCoord[4];

	float xmin, xmax, ymin, ymax, zmin, zmax;

    list<float> x1;
    list<float> x2;
    list<float> y1;
    list<float> y2;
    list<float> z1;
    list<float> z2;

    //* skip the first line
	getline(ifile, lineTemp);

	while( getline(ifile, lineTemp) ){

        if (lineTemp.empty()) continue;

		std::stringstream stringTokenizer(lineTemp);

        // Each panel line starts with 'Q' followed by 'C' (conductor) or 'D' (dielectric)
		stringTokenizer >> charTemp;
		if (charTemp != 'Q'){
            cerr << "ERROR: only 'Q' panel type is supported in this solver!" << endl;
            exit(1); 
		}

        stringTokenizer >> CondOrDiel;

        if(CondOrDiel == 'C'){
            this->isDielectric.push_back(false);
        }
        else if(CondOrDiel == 'D'){
            this->isDielectric.push_back(true);
            DielPanels++;
        }
        else{
            cout << "Invalid panel type '" << CondOrDiel << "'. Expected 'C' or 'D'.\n";
            exit(1);
        }

        // Group/conductor name
		stringTokenizer >> stringTemp;

        // Four corner vertices
		stringTokenizer >> xCoord[0] >> yCoord[0] >> zCoord[0]
					    >> xCoord[1] >> yCoord[1] >> zCoord[1]
					    >> xCoord[2] >> yCoord[2] >> zCoord[2]
					    >> xCoord[3] >> yCoord[3] >> zCoord[3];

        // Bounding box
		xmin = xmax = xCoord[0];
		ymin = ymax = yCoord[0];
		zmin = zmax = zCoord[0];
		for (int ci=1; ci<4; ci++){
			if ( xmin > xCoord[ci] ) xmin = xCoord[ci];
			if ( xmax < xCoord[ci] ) xmax = xCoord[ci];
			if ( ymin > yCoord[ci] ) ymin = yCoord[ci];
			if ( ymax < yCoord[ci] ) ymax = yCoord[ci];
			if ( zmin > zCoord[ci] ) zmin = zCoord[ci];
			if ( zmax < zCoord[ci] ) zmax = zCoord[ci];
		}
		x1.push_back(xmin); x2.push_back(xmax);
		y1.push_back(ymin); y2.push_back(ymax);
		z1.push_back(zmin); z2.push_back(zmax);

         // find the sign of the normal vector
        // Compute normal vector from vertices: (v1 - v0) x (v3 - v0)
        double ux = xCoord[1] - xCoord[0];
        double uy = yCoord[1] - yCoord[0];
        double uz = zCoord[1] - zCoord[0];
        double vx = xCoord[3] - xCoord[0];
        double vy = yCoord[3] - yCoord[0];
        double vz = zCoord[3] - zCoord[0];
        double nx = uy * vz - uz * vy;
        double ny = uz * vx - ux * vz;
        double nz = ux * vy - uy * vx;

        int nSign = 1;
        if (std::abs(nx) >= std::abs(ny) && std::abs(nx) >= std::abs(nz)) {
            nSign = (nx >= 0.0) ? 1 : -1;
        } else if (std::abs(ny) >= std::abs(nx) && std::abs(ny) >= std::abs(nz)) {
            nSign = (ny >= 0.0) ? 1 : -1;
        } else {
            nSign = (nz >= 0.0) ? 1 : -1;
        }
        this->normalSigns.push_back(nSign);

        if(CondOrDiel == 'D'){
            // New dielectric group?
            if(currentDielectricName.compare(stringTemp) != 0){ // check if we changed dielectric
                if(nDielPanelTemp > 0){
                    nDielPanelsList.push_back(nDielPanelTemp); //save the amount of the total panels of the previous dielectric group before resetting
                }
                currentDielectricName = stringTemp; //save the current dielectric group name 
                this->nDiels++; // amount of total dielectrics
                nDielPanelTemp = 0; // reset this counter to count the amount of panels in the current dielectric group
            }
            nDielPanelTemp++; // add to the panels of the current dielectric group 

            // Read relative permittivities: eps_a on the +normal side, eps_b on the -normal side
            // The normal direction is inferred from the panel orientation (panel facing direction).
            if(!(stringTokenizer >> epsA >> epsB)){
                cerr << "ERROR: Dielectric panel '" << currentDielectricName
                     << "' is missing permittivity values (eps_a eps_b)." << endl;
                exit(1);
            }
            this->materialA.push_back(epsA);
            this->materialB.push_back(epsB);
            this->nTotalDielPanels++;

        }
        else { // CondOrDiel == 'C' 
            // Conductor panels: permittivity = 1 (vacuum) on both sides
            this->materialA.push_back(1.0);
            this->materialB.push_back(1.0);

            if ( currentConductorName.compare(stringTemp) != 0 ){// check if we changed conductor
                // if(nWirePanelTemp>0){
                //     nWirePanelsList.push_back(nWirePanelTemp);
                // }
                currentConductorName = stringTemp; //save the current conductor group name 
                this->nWires++;// amount of total conudctors
                nWirePanelsList.push_back(nWirePanelTemp);//save the amount of the total panels of the previous conductor group before resetting
                nWirePanelTemp = 0;// reset this counter to count the amount of panels in the current conductor group
            }
            nWirePanelTemp++;// add to the panels of the current conductor group
            this->nTotalWirePanels++;
        }

        this->nPanels++; //count the total amount of panels conductor or dielectric
	}

    // After the while loop, register remaining groups
    if (nWirePanelTemp > 0) {
        nWirePanelsList.push_back(nWirePanelTemp);
    }
    if (nDielPanelTemp > 0) {
        nDielPanelsList.push_back(nDielPanelTemp);
    }



    // nCoefs = total unknowns (conductor charges + dielectric polarization charges)
	this->nCoefs = this->nTotalWirePanels;

    // Allocate and fill per-conductor panel count arrays
	this->nWirePanels = new int[this->nWires];
	this->nWireCoefs  = new int[this->nWires];

	int i = 0;
    // The list has a leading zero sentinel from the first push_back(0) before any conductor
	for ( std::list<int>::iterator it = ++nWirePanelsList.begin();
			it != nWirePanelsList.end(); ++it){
		this->nWirePanels[i] = *it;
		this->nWireCoefs [i] = *it;
		i++;
	}

    // Allocate and fill per-dielectric panel count arrays
    this->nDielPanels = new int[this->nDiels > 0 ? this->nDiels : 1];
    this->nDielCoefs  = new int[this->nDiels > 0 ? this->nDiels : 1];
    int j = 0;
    for (std::vector<int>::iterator it = nDielPanelsList.begin();
         it != nDielPanelsList.end(); ++it) {
        this->nDielPanels[j] = *it;
        this->nDielCoefs [j] = *it;
        j++;
    }

    cout << "Number of conductors        : " << this->nWires << endl;
    cout << "Number of dielectric groups : " << this->nDiels  << endl;


    //* this stores the exact geometry of each panel
    //indexes:
    //first: panel number
    //second: coordinate direction X,Y,Z
    //third: values MIN, MAX, LENGTH, CENTER
	this->panels 	= new float	[this->nPanels][3][4];

    //* this stores the orientation of each panel, the direction the panel is facing
	this->dirs 		= new int	[this->nPanels]; //NORMAL DIRECTION

    //* this stores the area of each panel
	this->areas		= new float	[this->nPanels];

	for ( int i=0; i<this->nPanels; i++){

        //* take the x,y,z values from the previous list, find the center in each axis and then store them in the panels[][][] array
		this->panels[i][X][MIN] = x1.front(); x1.pop_front();
		this->panels[i][X][MAX] = x2.front(); x2.pop_front();
		this->panels[i][X][LENGTH]
			   = this->panels[i][X][MAX] - this->panels[i][X][MIN];
		this->panels[i][X][CENTER]
			   = (this->panels[i][X][MAX] + this->panels[i][X][MIN])/2;

		this->panels[i][Y][MIN] = y1.front(); y1.pop_front();
		this->panels[i][Y][MAX] = y2.front(); y2.pop_front();
		this->panels[i][Y][LENGTH]
			   = this->panels[i][Y][MAX] - this->panels[i][Y][MIN];
		this->panels[i][Y][CENTER]
			   = (this->panels[i][Y][MAX] + this->panels[i][Y][MIN])/2;

		this->panels[i][Z][MIN] = z1.front(); z1.pop_front();
		this->panels[i][Z][MAX] = z2.front(); z2.pop_front();
		this->panels[i][Z][LENGTH]
			   = this->panels[i][Z][MAX] - this->panels[i][Z][MIN];
		this->panels[i][Z][CENTER]
			   = (this->panels[i][Z][MAX] + this->panels[i][Z][MIN])/2;



        //* find in which orientation the panel faces store it and then find each area and store it
        // find it by checking the length of each coordinate, example X-length=0, panel lies in the YZ plane, normal at X
		if ( std::abs(this->panels[i][X][LENGTH]) < zero ){ // in x-dir
			this->dirs[i] = X;
			this->areas[i] = this->panels[i][Y][LENGTH] * this->panels[i][Z][LENGTH];
		}else if ( std::abs(this->panels[i][Y][LENGTH]) < zero ){ // in y-dir
			this->dirs[i] = Y;
			this->areas[i] = this->panels[i][Z][LENGTH] * this->panels[i][X][LENGTH];
		}else{ // in z-dir
			this->dirs[i] = Z;
			this->areas[i] = this->panels[i][X][LENGTH] * this->panels[i][Y][LENGTH];
		}

	}

    //* For a uniform treatment of generateRHS
    //* indexIncrements describes how many panels belong to a basis function (1 new basis function, 0 the same basis function) for Arches
    //* here since the FasterCap format does not support Arches all the values of this array are one and all the panels are stored as FLAT ('F')
	this->indexIncrements 	= new int[this->nPanels];

    //*stores what kind of panel we have Fkat ('F') or Arch ('A')
	this->basisTypes		= new char[this->nPanels];

	for (int i=0; i<this->nPanels; i++){
		this->indexIncrements[i] = 1;
		this->basisTypes[i] = 'F';
	}


    //* Allocate dummy array for a uniform clear operation
	this->basisDirs			= new int[this->nPanels];
	for ( int i=0; i< this->nPanels; i++ ){
		this->basisDirs[i] = FLAT;
	}
	this->basisZs			= new float[1];
	this->basisShifts		= new float[1];

    //* Set flags
	this->isSolved = false;
	this->isLoaded = true;

	ifile.close();

    #ifdef DEBUG_ASPECT_RATIO_VALIDITY
    if( this->isPanelAspectRatioValid() == false ){
        cerr << "ERROR: Panel aspect ratio of FASTCAP file is not valid." << endl;
        exit(1);
    }
    #endif

    #ifdef ENABLE_DEBUG_PRINTS
    printAllData();
    // exit(2);
    #endif


}

#ifdef ENABLE_DEBUG_PRINTS
void Caplet::printAllData() const {
    using namespace std;

    cout << "\n========== CAPLET DATA DUMP ==========\n";
    cout << "isLoaded: " << isLoaded << endl;
    cout << "isSolved: " << isSolved << endl;
    cout << "nWires (conductors): " << nWires << endl;
    cout << "nDiels (dielectric groups): " << nDiels << endl;
    cout << "nPanels (total): " << nPanels << endl;
    cout << "nWirePanels (dielectric panels): " << nTotalWirePanels << endl;
    cout << "nDielPanels (dielectric panels): " << nTotalDielPanels << endl;
    cout << "nCoefs (basis functions, = conductor panels): " << nCoefs << endl;


    // Per-conductor panel counts
    cout << "\n--- Conductor groups ---" << endl;
    for (int i = 0; i < nWires; ++i) {
        cout << "Wire " << i << " : nWirePanels=" << nWirePanels[i]
             << ", nWireCoefs=" << nWireCoefs[i] << endl;
    }

    // Per-dielectric group counts (if any)
    if (nDiels > 0) {
        cout << "\n--- Dielectric groups ---" << endl;
        for (int i = 0; i < nDiels; ++i) {
            cout << "Dielectric group " << i << " : nDielPanels=" << nDielPanels[i]
                 << ", nDielCoefs=" << nDielCoefs[i] << endl;
        }
    }

    // Panel details
    cout << "\n--- Panel data (index, type, geometry, normal, area, eps1, eps2) ---" << endl;
    for (int i = 0; i < nPanels; ++i) {
        // Panel type
        const char* typeStr = (i < (int)isDielectric.size() && isDielectric[i]) ? "DIELECTRIC" : "CONDUCTOR";

        cout << "Panel " << i << " : " << typeStr << endl;
        cout << "  X: min=" << panels[i][X][MIN] << " max=" << panels[i][X][MAX]
             << " len=" << panels[i][X][LENGTH] << " center=" << panels[i][X][CENTER] << endl;
        cout << "  Y: min=" << panels[i][Y][MIN] << " max=" << panels[i][Y][MAX]
             << " len=" << panels[i][Y][LENGTH] << " center=" << panels[i][Y][CENTER] << endl;
        cout << "  Z: min=" << panels[i][Z][MIN] << " max=" << panels[i][Z][MAX]
             << " len=" << panels[i][Z][LENGTH] << " center=" << panels[i][Z][CENTER] << endl;

        // Normal direction
        const char* dirName = "?";
        if (dirs[i] == X) dirName = "X";
        else if (dirs[i] == Y) dirName = "Y";
        else if (dirs[i] == Z) dirName = "Z";
        cout << "  normal dir = " << dirName << " (value " << dirs[i] << ")" << endl;
        cout << "  area = " << areas[i] << endl;

        // Permittivities (materialA/B)
        if (i < (int)materialA.size() && i < (int)materialB.size()) {
            cout << "  eps1 (side of normal) = " << materialA[i] << endl;
            cout << "  eps2 (opposite side)  = " << materialB[i] << endl;
        } else {
            cout << "  WARNING: materialA/B not available for this panel" << endl;
        }

        // Basis function info (always flat for FastCap files)
        cout << "  basisType = " << basisTypes[i] << " (F=flat)" << endl;
        cout << "  indexIncrement = " << indexIncrements[i] << endl;
        cout << "  basisDir = " << basisDirs[i] << " (FLAT=" << FLAT << ")" << endl;
        cout << "  basisZ = " << basisZs[i] << endl;
        cout << "  basisShift = " << basisShifts[i] << endl;

        cout << "  ---" << endl;
    }

    cout << "========================================\n" << endl;
}
#endif


void Caplet::loadCapletFile(const std::string filename){
    //* Caplet format
	//* Line 1: nWire
	//* Line 2: nShape for each wire
	//* Line 3: nTotalShape
    //* Line n: shape description lines, 12 numbers per line
	//* shapeType, indexIncrement, XL, XU, YL, YU, ZL, ZU, dir, shapeDir, shapeNormalDistance, shapeShift
	//* shapeType: either F for Flat, A for Arch, and S for Side arch
	//* indexIncrement: 
	//*   a basis function can consists of multiple shapes. The value here can be 0 or 1.
	//*   1: the beginning of a new basis function
	//*   0: not increment, meaning the current shape is combined with the shape before with value 1
	//* XL, XU, YL, YU, ZL, ZU: the range limit of each direction
	//* dir: normal direction of the rectangle that supports the shape
	//*   0 for x-dir, 1 for y-dir, 2 for z-dir
	//* shapeDir: shape varying direction
	//*   0 for x-dir, 1 for y-dir, 2 for z-dir
	//* shapeNormalDistrance: signed distance between the support rectangle and the neighborhoold rectangle
	//*   positive: decaying in the positive direction
	//*   negative: decaying in the negative direction
	//* shapeShift: shift parameter that indicates the shift of starting point of 



	this->clear();
    this->isInstantiable = true;

    ifstream ifile(filename.c_str());
    //* Check if open successfully
	if (!ifile){
        cerr << "ERROR: cannot open the file: " << filename << endl;
		return;
	}

	ifile >> this->nWires;
	this->nWirePanels 		= new int[this->nWires];
	this->Cmat				= new float[this->nWires*this->nWires];
	for ( int i=0; i<this->nWires; i++){
		ifile >> this->nWirePanels[i];
	}
	ifile >> this->nPanels;
	this->panels			= new float	[this->nPanels][3][4];
	this->dirs				= new int	[this->nPanels];
	this->areas	 			= new float	[this->nPanels];
	this->indexIncrements 	= new int	[this->nPanels];
	this->basisTypes		= new char	[this->nPanels];
	this->basisDirs			= new int	[this->nPanels];
	this->basisZs			= new float	[this->nPanels];
	this->basisShifts		= new float	[this->nPanels];

	this->nCoefs = 0;

	for (int i=0; i<this->nPanels; i++){
		ifile >> this->basisTypes[i];
		ifile >> this->indexIncrements[i];
		this->nCoefs += this->indexIncrements[i];

		ifile >> this->panels[i][X][MIN];
		ifile >> this->panels[i][X][MAX];
		this->panels[i][X][LENGTH]
			= this->panels[i][X][MAX] - this->panels[i][X][MIN];
		this->panels[i][X][CENTER]
			= (this->panels[i][X][MAX]+this->panels[i][X][MIN])/2;

		ifile >> this->panels[i][Y][MIN];
		ifile >> this->panels[i][Y][MAX];
		this->panels[i][Y][LENGTH]
			= this->panels[i][Y][MAX] - this->panels[i][Y][MIN];
		this->panels[i][Y][CENTER]
			= (this->panels[i][Y][MAX]+this->panels[i][Y][MIN])/2;

		ifile >> this->panels[i][Z][MIN];
		ifile >> this->panels[i][Z][MAX];
		this->panels[i][Z][LENGTH]
			= this->panels[i][Z][MAX] - this->panels[i][Z][MIN];
		this->panels[i][Z][CENTER]
			= (this->panels[i][Z][MAX]+this->panels[i][Z][MIN])/2;

		ifile >> this->dirs[i];
		switch (this->dirs[i]){
		case X:
			this->areas[i] = this->panels[i][Y][LENGTH] * this->panels[i][Z][LENGTH];
			break;
		case Y:
			this->areas[i] = this->panels[i][Z][LENGTH] * this->panels[i][X][LENGTH];
			break;
		case Z:
			this->areas[i] = this->panels[i][X][LENGTH] * this->panels[i][Y][LENGTH];
		}

		ifile >> this->basisDirs[i];
		ifile >> this->basisZs[i];
		ifile >> this->basisShifts[i];
	}


	// this->P 	= new float[this->nCoefs*this->nCoefs];
	// this->rhs 	= new float[this->nCoefs*this->nWires];
	// this->coefs = new float[this->nCoefs*this->nWires];

	// construct nWireCoefs
	this->nWireCoefs = new int[this->nWires];
	int ind = -1;
	for (int i=0; i<this->nWires; i++){
		this->nWireCoefs[i] = 0;
		for (int j=0; j<this->nWirePanels[i]; j++){
			ind++;
			this->nWireCoefs[i] += this->indexIncrements[ind];
		}
	}

	this->isSolved = false;
	this->isLoaded = true;

	ifile.close();

    #ifdef DEBUG_ASPECT_RATIO_VALIDITY
    if( this->isPanelAspectRatioValid() == false ){
        cerr << "ERROR: Panel aspect ratio of CAPLET file is not valid." << endl;
        exit(1);
    }
    #endif
}


void Caplet::saveCmat(const std::string filename){
	if (this->isSolved){
		std::ofstream ofile(filename.c_str());
		if (!ofile.is_open()){
            cerr << "ERROR: cannot write Cmat file: " << filename << endl;
			return;
		}
		//* DOUBLE block was commented before. Double check this.
        switch(this->mode){
        case FAST_GALERKIN:
            for (int i=0; i<this->nWires; i++){
                for (int j=0; j<this->nWires; j++){
                    ofile << this->Cmat[i + this->nWires * j] << " ";
                }
                ofile << endl;
            }
            break;
		case DOUBLE_GALERKIN:            
        case DOUBLE_COLLOCATION:
            for (int i=0; i<this->nWires; i++){
                for (int j=0; j<this->nWires; j++){
                    ofile << this->dCmat[i + this->nWires * j] << " ";
                }
                ofile << endl;
            }
            break;
        default:
            break;
        }

		ofile.close();
	}
}


void Caplet::saveCoefs(const std::string filename){
	if (this->isSolved){
        ofstream ofile(filename.c_str());
		if (!ofile){
            cerr << "ERROR: cannot write the file: " << filename << endl;
			return;
		}
        switch(this->mode){
        case FAST_GALERKIN:
            for (int i=0; i<this->nCoefs; i++){
                for (int j=0; j<this->nWires; j++){
                    ofile << this->coefs[i + this->nCoefs * j]  * 4*pi*epsilon0 << " ";
                }
                ofile << endl;
            }
            break;
        case DOUBLE_GALERKIN:
        case DOUBLE_COLLOCATION:
            for (int i=0; i<this->nCoefs; i++){
                for (int j=0; j<this->nWires; j++){
                    ofile << this->dcoefs[i + this->nCoefs * j] * 4*pi*epsilon0 << " ";
                }
                ofile << endl;
            }
            break;
        default:
            ;
        }

		ofile.close();
	}
}

//function to implement the dielectric interactions for the P matrix
//Pij(diel)=(Aj(ri-rj) . ni)/4*pi*epsilon0*dij^3
double Caplet::calCollocationDielectricEntry(int panel_i, int panel_j) {
    const double four_pi_eps0 = 4.0*M_PI*(double)epsilon0;

    // Self term (i == j): -lambda_i = -(eps_in + eps_out) / (2 * eps0 * (eps_in - eps_out))
    if (panel_i==panel_j) {
        double eps_in  = this->materialA[panel_i];
        double eps_out = this->materialB[panel_i];
        return -((eps_in + eps_out) / (2.0 * (double)epsilon0 * (eps_in - eps_out)));
    }

    // Off-diagonal (i != j): 1-point far-field approximation
    // Pij^diel = A_j * ((r_i - c_j) . n_hat_i) / (4 * pi * eps0 * d_ij^3)
    double dx = panels[panel_i][X][CENTER] - panels[panel_j][X][CENTER];
    double dy = panels[panel_i][Y][CENTER] - panels[panel_j][Y][CENTER];
    double dz = panels[panel_i][Z][CENTER] - panels[panel_j][Z][CENTER];
    double dij = std::sqrt(dx * dx + dy * dy + dz * dz);

    if (dij < 1e-18) 
        return 0.0;

    int ni_axis = dirs[panel_i];
    double ni_sign = (double)normalSigns[panel_i];

    double r_dot_n = 0.0;
    if (ni_axis == X)      
        r_dot_n=dx*ni_sign;
    else if (ni_axis == Y)
        r_dot_n=dy*ni_sign;
    else                   
        r_dot_n=dz*ni_sign;

    double Aj = (double)areas[panel_j];
    return (Aj*r_dot_n) / (four_pi_eps0*dij*dij*dij);
}

void Caplet::extractC(MODE mode){
	this->mode = mode;

    #ifdef CAPLET_TIMER
	this->fillingTime = 0;
	this->solvingTime = 0;
	this->totalTime = 0;
    #endif

    #ifdef CAPLET_INIT_ATAN_LOG
    caplet::atan(1);
    caplet::log(1);
    #endif

    //* Init MPI
    // MPI::Init();
    // int rank = MPI::COMM_WORLD.Get_rank();


    MPI_Init(nullptr, nullptr);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //* Subdivide panels if aspect ratio is too large
    this->modifyPanelAspectRatio();

    if ( rank==0 ){
		std::cout << "Number of conductors        : " << this->nWires << std::endl;
		std::cout << "Number of basis functions   : " << this->nCoefs << std::endl;
        std::cout << "Number of basis shapes      : " << this->nPanels << std::endl;
    }
	for (int iter = 0; iter < N_ITER; iter++){
        switch ( mode ){
    	case DOUBLE_GALERKIN:
    		this->extractCGalerkinDouble();
    		break;
        case FAST_GALERKIN:
			this->extractCGalerkin();
			break;
        case DOUBLE_COLLOCATION:
            this->extractCCollocationDouble();
            break;
		}
	}
    this->isSolved = true;

    // MPI::Finalize();
    MPI_Finalize();
    if( rank!= 0 ){
		return;
	}


    #ifdef CAPLET_TIMER
    std::cout << "Total extraction time (s)   : " << (this->totalTime)/N_ITER << std::endl;
    std::cout << "  Setup time (s)            : " << (this->fillingTime)/N_ITER << std::endl;
    std::cout << "  Solving time (s)          : " << (this->solvingTime)/N_ITER << std::endl;
    #endif

    //* Print Cmat
    std::cout << "Cmat" << std::endl;
    this->printCmat();
	std::cout << endl;

	#ifdef DEBUG_PRINT_P
    std::cout << "Store system matrix in 'pmatrix' and rhs in 'rhs'" << std::endl;
    ofstream pmatrixOut("pmatrix");
    printP(pmatrixOut);
    pmatrixOut.close();
    ofstream rhsOut("rhs");
    printRHS(rhsOut);
    rhsOut.close();
    #endif
}



//***************************
//*
//* DOUBLE COOLLOCATION MODE
//*
//*
void Caplet::extractCCollocationDouble(){

    //find the total amount of groups
    int nTotalGroups = this->nWires + this->nDiels;
    int nTotalUnknowns = nCoefs + nTotalDielPanels;

    int rank = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if ( this->isLoaded == true ){
        //* Allocate system memory for double precision
        if (rank==0){
            this->dP 	= new double[nTotalUnknowns*nTotalUnknowns];
            this->drhs 	= new double[nTotalUnknowns*this->nWires];
            this->dcoefs= new double[nTotalUnknowns*this->nWires];
            this->dCmat = new double[this->nWires*this->nWires];
        }
        else{
            this->dP 	= new double[1];
            this->drhs 	= new double[1];
            this->dcoefs= new double[1];
            this->dCmat = new double[1];
            return;
        }

	}else{
		std::cerr << "ERROR: structure file is not yet loaded" << std::endl;
	}
    #ifdef CAPLET_TIMER
    this->timeStart = MPI_Wtime();
    #endif

    this->generateCollocationPMatrixDouble();
    this->generateRHSDouble();

    #ifdef CAPLET_TIMER
    this->timeAfterFilling = MPI_Wtime();
    #endif


    //* Solve the system
	int* 	ipiv = new int[nTotalUnknowns];
	int  	info;

    // int info;
    dgesv_(&nTotalUnknowns, &this->nWires, dP, &nTotalUnknowns, ipiv,
        dcoefs, &nTotalUnknowns, &info);
    if (info != 0) {
        std::cout << "ERROR: dgesv returned info = " << info << std::endl;
    }

	delete[] ipiv;

	char 	transA 	= 't';
	char 	transB 	= 'n';
	double 	alpha 	= 4*pi*epsilon0;
	double 	beta 	= 0.0;

    //* Use matrix-matrix product to compute Cmat from coefs
    dgemm_(&transA, &transB,
            &this->nWires, &this->nWires, &this->nCoefs,
            &alpha, this->drhs, &this->nCoefs,
            this->dcoefs, &this->nCoefs,
            &beta, this->dCmat, &this->nWires);

    #ifdef CAPLET_TIMER
    this->timeAfterSolving = MPI_Wtime();
 	this->fillingTime 	+= this->timeAfterFilling - this->timeStart;
 	this->solvingTime 	+= this->timeAfterSolving - this->timeAfterFilling;
 	this->totalTime		+= this->timeAfterSolving - this->timeStart;
    #endif
}


void Caplet::generateCollocationPMatrixDouble() {

    // Zero out dP
    double zero = 0.0;
    int    inc  = 1;
    int    nC   = nCoefs * nCoefs;
    dscal_(&nC, &zero, dP, &inc);
    double Q;

    //if there are no dielectric panels just use the the original code for conductors
    if(this->nTotalDielPanels==0){
        // Build index vector: ind[i] maps panel index i to its row/col in dP.
        // (For FLAT panels, indexIncrements[i] = 1 always, so ind[i] = i.)
        int* ind = new int[nPanels];
        ind[0] = 0;
        for (int i = 1; i < nPanels; i++) {
            ind[i] = ind[i-1] + indexIncrements[i];
        }

        const int nK = nPanels * nPanels;

        #ifdef CAPLET_OPENMP
            #pragma omp parallel for num_threads(CAPLET_OPENMP_NUM_THREADS)
        #endif
        for (int k = 0; k < nK; k++) {
            int j = k / nPanels;   // source panel  (column)
            int i = k - j * nPanels; // test panel  (row)

            double entry = 0.0;

            // if ( isDielectric[i] ) {
            //     // --- Dielectric row: E-kernel ----------------------------------
            //     entry = calCollocationEEntryDouble(i, j);
            // } else {
                // --- Conductor row: P-kernel -----------------------------------
                // Original collocation formula (1/r potential integral).
            entry = calGalerkinPEntry(i, j);
            // }

            dP[ ind[i] + nCoefs * ind[j] ] += entry;
        }

        delete[] ind;  
        // cout << "HI" <<endl;
        #ifdef ENABLE_P_MATRIX_PRINT
            print_dP(this->nPanels);
        #endif

        return;
    }
    else{
        // extended case with dielectrics
        int nTotal = this->nCoefs + this->nTotalDielPanels;
        
        // Build mapping: panel index -> unknown index
        // so the panels are in the correct order for the matrix however they were read
        //create a vector with unitialised elements
        std::vector<int> panelToUnknown(this->nPanels, -1);
        int condIdx = 0, dielIdx = this->nTotalWirePanels;
        for (int i = 0; i < this->nPanels; ++i) {
            if (!isDielectric[i])
                panelToUnknown[i] = condIdx++;
            else
                panelToUnknown[i] = dielIdx++;
        }

        #ifdef ENABLE_DEBUG_PRINTS
        std::cout << "Panel to unknown mapping:" << std::endl;
        #endif
        for (int i = 0; i < nPanels; ++i) {
            #ifdef ENABLE_DEBUG_PRINTS
            std::cout << "Panel " << i << " (" << (isDielectric[i] ? "D" : "C") 
                    << ") -> unknown " << panelToUnknown[i] << std::endl;
            #endif
        }
            // Fill matrix
        // for (int eqPanel = 0; eqPanel < this->nPanels; ++eqPanel) {
        //      //find to which uknown the current panel corresponds to 
        //      //for example eqPanel=0 could correspond to dielectric which is after the conductors in the P matrix
        //     int row = panelToUnknown[eqPanel];
        //     if (row == -1)
        //         continue;

        //     if (!isDielectric[eqPanel]) {  // Conductor equation
        //         for (int srcPanel = 0; srcPanel < this->nPanels; ++srcPanel) {
        //             int col = panelToUnknown[srcPanel];
        //             if (col == -1) continue;
        //             double P = calGalerkinPEntry(eqPanel, srcPanel);
        //             dP[row + nTotal * col] += P;
        //             #ifdef ENABLE_DEBUG_PRINTS
        //                 std::cout << std::setw(14) << std::scientific << P << " " << endl;
        //                 if (std::isnan(P) || std::isinf(P)) {
        //                     std::cout << "NaN/Inf in P for eqPanel=" << eqPanel << " srcPanel=" << srcPanel << std::endl;
        //                 }
        //             #endif
        //         }
        //     } else {  // Dielectric equation
        //         #ifdef ENABLE_DEBUG_PRINTS
        //             cout << "materialA: " << this->materialA[eqPanel] << "materialB: " << this->materialB[eqPanel];
        //         #endif

        //         double deps = this->materialA[eqPanel] - this->materialB[eqPanel];
        //         for (int srcPanel = 0; srcPanel < this->nPanels; ++srcPanel) {
        //             if(eqPanel==srcPanel){
        //                 // Q=(deps)/(2*epsilon0);
        //                 Q=0;
        //             }
        //             else{
        //             Q = calCollocationQEntry(eqPanel, srcPanel);

        //             }
        //             int col = panelToUnknown[srcPanel];
        //             if (col == -1) continue;
        //             // cout << "HIII\n";
        //             dP[row + nTotal * col] += deps * Q;
        //             #ifdef ENABLE_DEBUG_PRINTS
        //                 std::cout << std::setw(14) << std::scientific << Q << " " << endl;
        //                 if (std::isnan(Q) || std::isinf(Q)) {
        //                     std::cout << "NaN/Inf in P for eqPanel=" << eqPanel << " srcPanel=" << srcPanel << std::endl;
        //                 }
        //             #endif


        //             // cout << eqPanel, srcPanel, Q << endl;

        //         }
        //     }

        // }

        const double four_pi_eps0 = 4.0 * M_PI * (double)epsilon0;      

        for (int eqPanel = 0; eqPanel < this->nPanels; ++eqPanel) {
                int row = panelToUnknown[eqPanel];
                if (row == -1) continue;

                for (int srcPanel = 0; srcPanel < this->nPanels; ++srcPanel) {
                    int col = panelToUnknown[srcPanel];
                    if (col == -1) continue;

                    if (!isDielectric[eqPanel]) {
                        // Conductor row: Use the original calGalerkinPEntry as requested
                        double P_val = calGalerkinPEntry(eqPanel, srcPanel);
                        
                        // We divide by 4*pi*eps0 here so the conductor block's physical scale 
                        // perfectly matches the dielectric block's physical scale.
                        dP[ψολ + nTotal * col] += P_val / four_pi_eps0; 
                    } else {
                        // Dielectric row: Use the dielectric boundary condition function
                        dP[row + nTotal * col] += calCollocationDielectricEntry(eqPanel, srcPanel);
                    }
                }
            }

        #ifdef ENABLE_P_MATRIX_PRINT
            print_dP(nTotal, "dP_matrix_output.txt");
        #endif

    }
}

#ifdef ENABLE_P_MATRIX_PRINT
//PRINT IN TERMINAL
// void Caplet::print_dP(int nTotal) const {
//     std::cout << "\n=== Extended matrix dP (" << nTotal << "x" << nTotal << ") ===" << std::endl;
//     for (int i = 0; i < nTotal; ++i) { 
//         for (int j = 0; j < nTotal; ++j) {
//             // dP is column‑major: element (i,j) is stored at index  (1/(4*3.1415926535f*epsilon0)) i + nTotal*j
//             std::cout << std::setw(14) << std::scientific << dP[i + nTotal * j] << " ";
//         }
//         std::cout << std::endl;
//     }
//     std::cout << std::endl;
// }

//PRINT IN FILE
void Caplet::print_dP(int nTotal, const std::string& filename) const {
    // Open the file for writing (use std::ios::app if you want to append instead of overwrite)
    std::ofstream outFile(filename);
    
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    outFile << "\n=== Extended matrix dP (" << nTotal << "x" << nTotal << ") ===" << std::endl;
    for (int i = 0; i < nTotal; ++i) { 
        for (int j = 0; j < nTotal; ++j) {
            // dP is column-major: element (i,j) is stored at index i + nTotal*j
            outFile << std::setw(14) << std::scientific << dP[i + nTotal * j] << " ";
        }
        outFile << std::endl;
    }
    outFile << std::endl;
    
    outFile.close();
}
#endif

////RUUUNS
// double Caplet::calCollocationQEntry(int panel_i, int panel_j) {
//     // Get the rotated coordinates of the centroids
//     // We need to rotate panel_i so its normal becomes +Z,
//     // and apply the same rotation to panel_j.

//     // Temporary coordinate arrays for the two panels (like in calCollocationPEntryDouble)
//     float *coord_i[3][4];
//     float *coord_j[3][4];

//     // Choose rotation based on the normal direction of panel_i
//     //change locally the panels so their normals are at Z
//     switch (dirs[panel_i]) {
//         case X:
//             rotateX2Z(coord_i, panel_i);
//             rotateX2Z(coord_j, panel_j);
//             break;
//         case Y:
//             rotateY2Z(coord_i, panel_i);
//             rotateY2Z(coord_j, panel_j);
//             break;
//         case Z:
//             rotateZ2Z(coord_i, panel_i);
//             rotateZ2Z(coord_j, panel_j);
//             break;
//         default:
//             return 0.0;
//     }

//     cout << "panel_i: " << panel_i << " direction: " << dirs[panel_i];
//     cout << "panel_j: " << panel_j << " direction: " << dirs[panel_j];


//     // Now the centroids are in the rotated frame.
//     // For panel_i, its normal is +Z, so its centroid has coordinates (x_i, y_i, z_i).
//     // For panel_j, its centroid is (x_j, y_j, z_j) after the same rotation.
//     double xi = (*coord_i[X])[CENTER];
//     double yi = (*coord_i[Y])[CENTER];
//     double zi = (*coord_i[Z])[CENTER];
//     double xj = (*coord_j[X])[CENTER];
//     double yj = (*coord_j[Y])[CENTER];
//     double zj = (*coord_j[Z])[CENTER];

//     //compute the formula
//     double dx = xi - xj;
//     double dy = yi - yj;
//     double dz = zi - zj;
//     double r = sqrt(dx*dx + dy*dy + dz*dz);
//     double r3 = r * r * r;

//     // Avoid division by zero for self‑panel (i == j)
//     if (r3 < 1e-18) return 0.0;

//     // Q = A_j / (4π ε0) * ( - dz / r³ )
//     double Q = - areas[panel_j] * (dz / r3);

//     return Q;
// }

double Caplet::calCollocationQEntry(int panel_i, int panel_j) {
    // Vector from source panel j center to test panel i center
    double dx = panels[panel_i][X][CENTER] - panels[panel_j][X][CENTER];
    double dy = panels[panel_i][Y][CENTER] - panels[panel_j][Y][CENTER];
    double dz = panels[panel_i][Z][CENTER] - panels[panel_j][Z][CENTER];

    #ifdef ENABLE_DEBUG_PRINTS
        cout << "panel_i: " << panel_i << " direction: " << dirs[panel_i] << endl;
        cout << "panel_j: " << panel_j << " direction: " << dirs[panel_j] << endl;
    #endif

    // Component along the outward normal of panel i
    double ndotr = 0.0;
    switch (dirs[panel_i]) {
        case X: ndotr = dx; break;
        case Y: ndotr = dy; break;
        case Z: ndotr = dz; break;
        default: return 0.0;
    }

    double r = sqrt(dx*dx + dy*dy + dz*dz);
    double r3 = r * r * r;
    if (r3 < 1e-18) return 0.0;

    // Q = - (area of source) * (normal component) / r^3
    return - areas[panel_j] * ndotr / r3;
}

double Caplet::calCollocationEEntryDouble(int panel_d, int panel_j) {

    float *cd[3][4];   // coordinate pointers for dielectric panel d (test point)
    float *cj[3][4];   // coordinate pointers for source panel j

    const double eps_a = materialA[panel_d];  // permittivity on n_hat_d side
    const double eps_b = materialB[panel_d];  // permittivity on opposite side

    // --- Analytically exact: same-direction panels -----------------------
    // When panel_d and panel_j face the same direction, one rotation maps
    // both panels so that n_hat_d = Z.  calEColD then integrates Δz/r^3 over
    // panel j, where Δz is the displacement along n_hat_d.

    // if ( dirs[panel_d] == dirs[panel_j] ) {
    //     switch ( dirs[panel_d] ) {
    //     case X:
    //         this->rotateX2Z(cd, panel_d);
    //         this->rotateX2Z(cj, panel_j);
    //         break;
    //     case Y:
    //         this->rotateY2Z(cd, panel_d);
    //         this->rotateY2Z(cj, panel_j);
    //         break;
    //     default: // Z
    //         this->rotateZ2Z(cd, panel_d);
    //         this->rotateZ2Z(cj, panel_j);
    //         break;
    //     }
    //     return -(eps_a - eps_b) * calEColD(cd, cj);


    // }

    cout << "panel_d: " << panel_d << "direction: " << dirs[panel_d];
    cout << "panel_j: " << panel_j << "direction: " << dirs[panel_j];


                // --- Far-field dipole: cross-direction panels --------------------------
    // For panels with different orientations (e.g., vertical sidewall + horizontal
    // dielectric layer), the outward normal of d does not align with the natural
    // integration direction of j.  We use the monopole/dipole approximation:
    //
    //   d/dn_hat_d (1/r)  ≈  -(n_hat_d · (r_d - r_j)) / |r_d - r_j|^3
    //
    // accurate when |r_d - r_j| >> max(size_d, size_j).
 
    const double dx = panels[panel_d][X][CENTER] - panels[panel_j][X][CENTER];
    const double dy = panels[panel_d][Y][CENTER] - panels[panel_j][Y][CENTER];
    const double dz = panels[panel_d][Z][CENTER] - panels[panel_j][Z][CENTER];
    const double r2 = dx*dx + dy*dy + dz*dz;
    if (r2 < 1e-30) return 0.0;   // coincident panel centers (guard)
    const double r3 = r2 * std::sqrt(r2);
 
    // Normal component of (r_d - r_j) along n_hat_d
    double ndotr = 0.0;
    switch ( dirs[panel_d] ) {
    case X: ndotr = dx; break;
    case Y: ndotr = dy; break;
    default: ndotr = dz; break;
    }
 
    // E_dj ≈ -(eps_a - eps_b) * (-ndotr / r^3) = (eps_a - eps_b) * ndotr / r^3
    // (the leading minus comes from d/dn(1/r) = -(n · (r_d - r'))/r^3,
    //  here evaluated at r' = r_j center)
    return -(eps_a - eps_b) * (-ndotr / r3);


}


double Caplet::calCollocationPEntryDouble(int panel1, int panel2){
    //* Test point : panel1.center
    //  Integration: panel2

    //* A nDim-element array of pointers pointing to a nBit-element array
    float *coord_ptr_1[3][4];
    float *coord_ptr_2[3][4];

    //* Rotate, mirror, and call proper integrals
    switch( dirs[panel2] ){
    case X:
        this->rotateX2Z(coord_ptr_1, panel1);
        this->rotateX2Z(coord_ptr_2, panel2);
        return calColD( coord_ptr_1, coord_ptr_2);
        break;
    case Y:
        this->rotateY2Z(coord_ptr_1, panel1);
        this->rotateY2Z(coord_ptr_2, panel2);
        return calColD( coord_ptr_1, coord_ptr_2);
        break;
    case Z:
        return calColD( coord_ptr_1, coord_ptr_2);
        break;
    }
    return 0.0; // dummy return
}


void Caplet::generateRHSDouble() {


    if (nTotalDielPanels == 0) {
            int inc  = 1;
        int nRHS = this->nWires * this->nCoefs;
        // Zero out dcoefs — it will hold the RHS initially
        double alpha = 0.0;
        dscal_(&nRHS, &alpha, this->dcoefs, &inc);

        int rowIndex    = -1;   // current row in dP (incremented by indexIncrements)
        int columnIndex =  0;   // current wire index
        int coefCounter =  0;   // counts coefs assigned to current wire

        for (int i = 0; i < this->nPanels; i++) {

            // Always advance the row index (even for dielectric panels,
            // since they occupy rows in the matrix)
            rowIndex += this->indexIncrements[i];

            // Skip the conductor coef-counter / column-index logic entirely.
            if ( isDielectric[i] ) {
                continue;
            }

            // --- Conductor panels: assign area to the correct wire column ------
            coefCounter += this->indexIncrements[i];

            // Advance to next wire once we have filled all its coefs
            if ( coefCounter > this->nWireCoefs[columnIndex] ) {
                coefCounter = 1;
                columnIndex++;
            }

            // Compute the effective "area" for this basis function.
            // For FLAT panels this is just areas[i].
            // For ARCH/SIDE panels the shape-weighted area is computed below.
            float areaVal = this->areas[i];
            float (*func)(float x, float w) = &arch;
            switch (this->basisTypes[i]) {
            case 'S':
                func = &side;
                // FALL THROUGH to 'A' case intentionally
            case 'A': {
                float b   = this->panels[i][this->basisDirs[i]][LENGTH];
                float xm  = b / 2;
                float xr  = b / 2;
                int dir = (this->basisDirs[i] + 1) % 3;
                if (dir == dirs[i]) dir = (dir + 1) % 3;
                float w   = this->panels[i][dir][LENGTH];

                float val = 0;
                int init_i = 0;
                if ((gauss_n % 2) == 1) {
                    init_i = 1;
                    val += (*gauss::w[gauss_n])[0] * func(xm, w);
                }
                int gauss_n2 = (gauss_n + 1) / 2;
                for (int gi = init_i; gi < gauss_n2; gi++) {
                    float dx = xr * (*gauss::p[gauss_n])[gi];
                    val += (*gauss::w[gauss_n])[gi]
                        * (func(xm + dx, w) + func(xm - dx, w));
                }
                areaVal = val * xr * w;
                this->areas[i] = areaVal;   // update cached area
                break;
            }
            default:
                break;
            }

            this->dcoefs[ columnIndex * this->nCoefs + rowIndex ] += areaVal;
        }

    // Copy RHS into drhs for safekeeping (dcoefs will be overwritten by solver)
        dcopy_(&nRHS, this->dcoefs, &inc, this->drhs, &inc);
    }

     // ----- With dielectrics: build extended RHS -----
    int nTotal = this->nCoefs + this->nTotalDielPanels;

             int inc  = 1;
    int nRHS = this->nWires * nTotal;
    double zero = 0.0;
    dscal_(&nRHS, &zero, this->dcoefs, &inc);


    #ifdef ENABLE_DEBUG_PRINTS
        std::cout << "DEBUG: nTotal=" << nTotal << ", nWires=" << nWires << ", nRHS=" << nRHS << std::endl;
        std::cout << "DEBUG: nCoefs=" << nCoefs << ", nTotalDielPanels=" << nTotalDielPanels << std::endl;
        std::cout << "DEBUG: nPanels=" << nPanels << std::endl;

        // Print nWirePanels
        int sumWire = 0;
        for (int i=0; i<nWires; ++i) {
            std::cout << "nWirePanels[" << i << "]=" << nWirePanels[i] << std::endl;
            sumWire += nWirePanels[i];
        }
        std::cout << "Sum nWirePanels = " << sumWire << ", nCoefs = " << nCoefs << std::endl;
    #endif


    

  // Build mapping for conductor panels only
    std::vector<int> condRow(nPanels, -1);      // row index in dcoefs for each conductor panel
    std::vector<int> condWire(nPanels, -1);     // wire index for each conductor panel

    int currentRow = 0;
    int currentWire = 0;
    int panelsSeenInWire = 0;
    for (int i = 0; i < nPanels; ++i) {
        if (isDielectric[i]) continue;          // skip dielectrics
        condRow[i] = currentRow++;
        condWire[i] = currentWire;
        panelsSeenInWire++;
        if (panelsSeenInWire == nWirePanels[currentWire]) {
            currentWire++;
            panelsSeenInWire = 0;
        }
    }

    // Fill dcoefs with RHS values (one column per excitation)
    for (int exc = 0; exc < nWires; ++exc) {
        for (int i = 0; i < nPanels; ++i) {
            if (isDielectric[i]) continue;
            int row = condRow[i];
            int wid = condWire[i];
            double val = (wid == exc) ? 1.0 : 0.0;
            // Correct column‑major indexing: index = row + nTotal * exc
            dcoefs[row + nTotal * exc] = val;
        }
    }

    // Copy to drhs (same layout)
    dcopy_(&nRHS, this->dcoefs, &inc, this->drhs, &inc);

    #ifdef ENABLE_DEBUG_PRINTS
    std::cout << "RHS matrix drhs:" << std::endl;
    for (int i = 0; i < nTotal; ++i) {
        for (int j = 0; j < nWires; ++j) {
            std::cout << std::setw(12) << std::scientific << drhs[i + nTotal * j] << " ";
        }
        std::cout << std::endl;
    }
    #endif
}


//__________________________________________________________
//*
//* DOUBLE GALERKIN MODE
//*
void Caplet::extractCGalerkinDouble(){

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if ( this->isLoaded == true ){
        if( rank==0 ){
            this->dP = new double[this->nCoefs*this->nCoefs];
        }else{
            this->dP = new double[1];
        }

        this->drhs 	 = new double[this->nCoefs*this->nWires];
        this->dcoefs = new double[this->nCoefs*this->nWires];
        this->dCmat  = new double[this->nWires*this->nWires];

    }else{
        std::cerr << "ERROR: structure file is not yet loaded" << std::endl;
    }
    #ifdef CAPLET_TIMER
    this->timeStart = MPI_Wtime();
    #endif


    #ifdef CAPLET_MPI
    this->generateGalerkinPMatrixDoubleMPI();
    #else
    this->generateGalerkinPMatrixDouble();
    #endif

    int rank2;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank2!=0){
        return;
    }
    //* End of core with non-zero rank

    this->generateRHSDouble();

    #ifdef CAPLET_TIMER
    this->timeAfterFilling = MPI_Wtime();
    #endif


    //* Solve the system
    int  	info;
    char 	uplo = 'u';

    //* Query optimal workspace size
    double*	work = new double[1];
    int		lwork = -1;
    int*	ipiv = new int[this->nCoefs];
    dsysv_(&uplo, &nCoefs, &nWires, dP, &nCoefs, ipiv, this->dcoefs, &nCoefs, work, &lwork, &info);
    lwork = work[0];
    delete[] work;

    //* Solve system using optimal work length
    work = new double[lwork];
    dsysv_(&uplo, &nCoefs, &nWires, dP, &nCoefs, ipiv, this->dcoefs, &nCoefs, work, &lwork, &info);
    delete[] ipiv;
    delete[] work;


    //* Use matrix-matrix product to compute Cmat from coefs
    char 	transA 	= 't';
    char 	transB 	= 'n';
    double 	alpha 	= 4*pi*epsilon0;
    double 	beta 	= 0.0;
    dgemm_(&transA, &transB,
            &this->nWires, &this->nWires, &this->nCoefs,
            &alpha, this->drhs, &this->nCoefs,
            this->dcoefs, &this->nCoefs,
            &beta, this->dCmat, &this->nWires);

    #ifdef CAPLET_TIMER
    this->timeAfterSolving = MPI_Wtime();

    this->fillingTime 	+= this->timeAfterFilling - this->timeStart;
    this->solvingTime 	+= this->timeAfterSolving - this->timeAfterFilling;
    this->totalTime		+= this->timeAfterSolving - this->timeStart;
    #endif

}

//__________________________________________________________
//*
//* FAST GALERKIN MODE
//*
void Caplet::extractCGalerkin(){

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    if ( this->isLoaded == true ){
        if( rank==0 ){
            this->P 	= new float[this->nCoefs*this->nCoefs];
        }else{
            this->P		= new float[1];
        }

        this->rhs 	= new float[this->nCoefs*this->nWires];
        this->coefs = new float[this->nCoefs*this->nWires];
        this->Cmat  = new float[this->nWires*this->nWires];

    }else{
        std::cerr << "ERROR: structure file is not yet loaded" << std::endl;
    }
    #ifdef CAPLET_TIMER
    this->timeStart = MPI_Wtime();
    #endif


    #ifdef CAPLET_MPI
    this->generateGalerkinPMatrixMPI();
    #endif

    #ifndef CAPLET_MPI
    this->generateGalerkinPMatrix();
    #endif

    int rank2;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank2!=0){
        return;
    }
    this->generateRHS();


    #ifdef CAPLET_TIMER
    this->timeAfterFilling = MPI_Wtime();
    #endif


    //* Solve the system
    int  	info;
    char 	uplo = 'u';

    //* Query optimal workspace size
    float*	work = new float[1];
    int		lwork = -1;
    int*	ipiv = new int[this->nCoefs];
    ssysv_(&uplo, &nCoefs, &nWires, P, &nCoefs, ipiv, coefs, &nCoefs, work, &lwork, &info);
    lwork = work[0];
    delete[] work;

    //* Solve system using optimal work length
    work = new float[lwork];
    ssysv_(&uplo, &nCoefs, &nWires, P, &nCoefs, ipiv, coefs, &nCoefs, work, &lwork, &info);
    delete[] ipiv;
    delete[] work;


    //* Use matrix-matrix product to compute Cmat from coefs
    char 	transA 	= 't';
    char 	transB 	= 'n';
    float 	alpha 	= 4*pi*epsilon0;
    float 	beta 	= 0.0;
    sgemm_(&transA, &transB,
            &this->nWires, &this->nWires, &this->nCoefs,
            &alpha, this->rhs, &this->nCoefs,
            this->coefs, &this->nCoefs,
            &beta, this->Cmat, &this->nWires);

    #ifdef CAPLET_TIMER
    this->timeAfterSolving = MPI_Wtime();

    this->fillingTime 	+= this->timeAfterFilling - this->timeStart;
    this->solvingTime 	+= this->timeAfterSolving - this->timeAfterFilling;
    this->totalTime		+= this->timeAfterSolving - this->timeStart;
    #endif
}


void Caplet::generateGalerkinPMatrix(){

    float zero = 0.0f;
    int   inc  = 1;
    int   nC   = nCoefs*nCoefs;
    sscal_(&nC, &zero, P, &inc);

    //* Construct ind_vec from indexIncrements
    int* ind = new int[nPanels];
    ind[0] = 0;
    for ( int i=1; i<nPanels; i++){
        ind[i] = ind[i-1] + indexIncrements[i];
    }
    const int nK = nPanels*(nPanels+1)/2;

    #ifdef CAPLET_OPENMP
        #pragma omp parallel for num_threads(CAPLET_OPENMP_NUM_THREADS)
    #endif
    for ( int k=0; k < nK; k++ ){
        int j = int((sqrt(double(1+8*k))-1)/2);
        int i = k - j*(j+1)/2;

        float result = calGalerkinPEntry(i,j);

        if ( (i!=j) && (ind[i]==ind[j]) ){
            P[ ind[i] + nCoefs*ind[j] ] += result*2;
        }else{
            P[ ind[i] + nCoefs*ind[j] ] += result;
        }
    }

    delete[] ind;

    //* Apply dielectric boundary condition (float Galerkin version)
    int* ind2 = new int[nPanels];
    ind2[0] = 0;
    for ( int i=1; i<nPanels; i++){
        ind2[i] = ind2[i-1] + indexIncrements[i];
    }
    for ( int i = 0; i < nPanels; i++ ){
        if ( !isDielectric[i] ) continue;
        float epsA  = (float)this->materialA[i];
        float epsB  = (float)this->materialB[i];
        float kappa = (epsA - epsB) / (epsA + epsB);
        int ci = ind2[i];
        float pSelf = P[ ci + nCoefs*ci ];
        for ( int j = 0; j < nCoefs; j++ ){
            P[ ci + nCoefs*j ] *= -kappa;
        }
        P[ ci + nCoefs*ci ] = pSelf * 0.5f;
    }
    delete[] ind2;
}

void Caplet::generateGalerkinPMatrixDouble(){

    double zero = 0.0;
    int    inc  = 1;
    int    nC   = nCoefs*nCoefs;
    dscal_(&nC, &zero, dP, &inc);

    //* Construct ind_vec from indexIncrements
    int* ind = new int[nPanels];
    ind[0] = 0;
    for ( int i=1; i<nPanels; i++){
        ind[i] = ind[i-1] + indexIncrements[i];
    }
    const int nK = nPanels*(nPanels+1)/2;

    #ifdef CAPLET_OPENMP
        #pragma omp parallel for num_threads(CAPLET_OPENMP_NUM_THREADS)
    #endif
    for ( int k=0; k < nK; k++ ){
        int j = int((sqrt(double(1+8*k))-1)/2);
        int i = k - j*(j+1)/2;

        double result = calGalerkinPEntry(i,j);

        if ( (i!=j) && (ind[i]==ind[j]) ){
            dP[ ind[i] + nCoefs*ind[j] ] += result*2;
        }else{
            dP[ ind[i] + nCoefs*ind[j] ] += result;
        }
    }

    delete[] ind;

    //* Apply dielectric boundary condition (double Galerkin version)
    int* ind2 = new int[nPanels];
    ind2[0] = 0;
    for ( int i=1; i<nPanels; i++){
        ind2[i] = ind2[i-1] + indexIncrements[i];
    }
    for ( int i = 0; i < nPanels; i++ ){
        if ( !isDielectric[i] ) continue;
        double epsA  = this->materialA[i];
        double epsB  = this->materialB[i];
        double kappa = (epsA - epsB) / (epsA + epsB);
        int ci = ind2[i];
        double pSelf = dP[ ci + nCoefs*ci ];
        for ( int j = 0; j < nCoefs; j++ ){
            dP[ ci + nCoefs*j ] *= -kappa;
        }
        dP[ ci + nCoefs*ci ] = pSelf * 0.5;
    }
    delete[] ind2;
}

//* Convert lower triangular matrix (column major) index k to subscript i,j
inline void ltind2sub(int k, int& i, int& j){
    j = int((std::sqrt(double(1+8*k))-1)/2);
    i = k - j*(j+1)/2;
}
void Caplet::generateGalerkinPMatrixMPI(){

    int rank, numproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numproc);

    int totalK = nPanels*(nPanels+1)/2;
    int nK = totalK/numproc;

    int* startK 	= new int[numproc];
    int* lastK 		= new int[numproc];
    int* startC		= new int[numproc];
    int* lastC		= new int[numproc];

    for( int i=0; i<numproc; i++ ){
        startK[i] = i*nK;
        lastK[i] = (i+1)*nK-1;
    }
    lastK[numproc-1] = totalK-1;

    //* Construct ind_vec from indexIncrements
    int* ind = new int[nPanels];
    ind[0] = 0;
    for ( int i=1; i<nPanels; i++){
        ind[i] = ind[i-1] + indexIncrements[i];
    }

    int maxNC = -1;
    for( int r=0; r<numproc; r++ ){
        int i, startj, lastj;
        ltind2sub(startK[r], i, startj);
        ltind2sub(lastK[r], i,  lastj );
        startC[r] = ind[startj];
        lastC[r]  = ind[lastj];

        int nC = lastC[r]-startC[r]+1;
        if ( maxNC < nC ){
            maxNC = nC;
        }
    }

    float* ptrP;
    float* tempP;
    if ( rank==0 ){
        //* init this->P
        float zero = 0.0f;
        int   inc  = 1;
        int   nP   = nCoefs*nCoefs;
        sscal_(&nP, &zero, P, &inc);
        //* init tempP to cover the largest number of columns among all proccesses.
        int   nTempP = maxNC*nCoefs;
        tempP= new float[nTempP];
        ptrP = this->P;
    }else{
        float zero 	= 0.0f;
        int   inc  	= 1;
        int   nTempP= ( lastC[rank]-startC[rank]+1 )*nCoefs;
        tempP= new float[ nTempP ];
        sscal_(&nTempP, &zero, tempP, &inc);
        ptrP = tempP;
    }

    //* For each k index
    for ( int k = startK[rank] ; k <= lastK[rank] ; k++ ){
        //* Convert k index to i,j subscript
        int i,j;
        ltind2sub(k, i, j);
        //* Compute the P entry

        float result = calGalerkinPEntry(i,j);

        #ifdef DEBUG_DETECT_NAN_INF_ENTRY
        #include <cmath>
        if (isnan(result) || isinf(result)){
        	cout << "Detect nan or inf: (" << i << "," << j << ") = " << result << endl;
        	cout << "i: " 	<< panels[i][X][MIN] << ", " << panels[i][X][MAX] << ", " 
        		 			<< panels[i][Y][MIN] << ", " << panels[i][Y][MAX] << ", " 
        		 			<< panels[i][Z][MIN] << ", " << panels[i][Z][MAX] << ". Dir: "
        		 			<< dirs[i] << ", " << basisDirs[i] << endl;
        	cout << "j: " 	<< panels[j][X][MIN] << ", " << panels[j][X][MAX] << ", " 
        					<< panels[j][Y][MIN] << ", " << panels[j][Y][MAX] << ", " 	
        					<< panels[j][Z][MIN] << ", " << panels[j][Z][MAX] << ". Dir: "
        					<< dirs[j] << ", " << basisDirs[j] << endl;
        }
        #endif


        //* Combine rows or columns in place if consecutive panels
        //  belong to the same basis function
        if ( (i!=j) && (ind[i]==ind[j]) ){
            ptrP[ ind[i] + nCoefs*(ind[j]-startC[rank]) ] += result*2;
        }else{
            ptrP[ ind[i] + nCoefs*(ind[j]-startC[rank]) ] += result;
        }
    }

    //* Combine sub-matrices to rank0 node
    if( rank==0 ){
        MPI_Status status;
        for ( int i=1; i<numproc; i++ ){
            int copylen = (lastC[i]-startC[i]+1)*nCoefs;
            MPI_Recv(tempP, copylen, MPI_FLOAT, i, 0, MPI_COMM_WORLD, &status);
            float alpha = 1.0f;
            int inc = 1;
            ptrP = P+(startC[i]*nCoefs);
            saxpy_(&copylen, &alpha, tempP, &inc, ptrP, &inc);
        }
    }else{
        int copylen = (lastC[rank]-startC[rank]+1)*nCoefs;
        MPI_Send(tempP, copylen, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
    }


    delete[] startK;
    delete[] lastK;
    delete[] startC;
    delete[] lastC;
    delete[] tempP;
    delete[] ind;

}
void Caplet::generateGalerkinPMatrixDoubleMPI(){

    int rank, numproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numproc);

    int totalK = nPanels*(nPanels+1)/2;
    int nK = totalK/numproc;

    int* startK 	= new int[numproc];
    int* lastK 		= new int[numproc];
    int* startC		= new int[numproc];
    int* lastC		= new int[numproc];

    for( int i=0; i<numproc; i++ ){
        startK[i] = i*nK;
        lastK[i] = (i+1)*nK-1;
    }
    lastK[numproc-1] = totalK-1;

    //* Construct ind_vec from indexIncrements
    int* ind = new int[nPanels];
    ind[0] = 0;
    for ( int i=1; i<nPanels; i++){
        ind[i] = ind[i-1] + indexIncrements[i];
    }

    int maxNC = -1;
    for( int r=0; r<numproc; r++ ){
        int i, startj, lastj;
        ltind2sub(startK[r], i, startj);
        ltind2sub(lastK[r], i,  lastj );
        startC[r] = ind[startj];
        lastC[r]  = ind[lastj];

        int nC = lastC[r]-startC[r]+1;
        if ( maxNC < nC ){
            maxNC = nC;
        }
    }

    double* ptrP;
    double* tempP;
    if ( rank==0 ){
        //* init this->dP
        double zero = 0.0;
        int    inc  = 1;
        int    nP   = nCoefs*nCoefs;
        dscal_(&nP, &zero, dP, &inc);
        //* init tempP to cover the largest number of columns among all proccesses.
        int   nTempP = maxNC*nCoefs;
        tempP= new double[nTempP];
        ptrP = this->dP;
    }else{
        double zero = 0.0;
        int    inc  	= 1;
        int    nTempP= ( lastC[rank]-startC[rank]+1 )*nCoefs;
        tempP= new double[ nTempP ];
        dscal_(&nTempP, &zero, tempP, &inc);
        ptrP = tempP;
    }

    //* For each k index
    for ( int k = startK[rank] ; k <= lastK[rank] ; k++ ){
        //* Convert k index to i,j subscript
        int i,j;
        ltind2sub(k, i, j);
        //* Compute the P entry

        double result = calGalerkinPEntry(i,j);

        #ifdef DEBUG_DETECT_NAN_INF_ENTRY
        #include <cmath>
        if (isnan(result) || isinf(result)){
        	cout << "Detect nan or inf: (" << i << "," << j << ") = " << result << endl;
        	cout << "i: " 	<< panels[i][X][MIN] << ", " << panels[i][X][MAX] << ", " 
        		 			<< panels[i][Y][MIN] << ", " << panels[i][Y][MAX] << ", " 
        		 			<< panels[i][Z][MIN] << ", " << panels[i][Z][MAX] << ". Dir: "
        		 			<< dirs[i] << ", " << basisDirs[i] << endl;
        	cout << "j: " 	<< panels[j][X][MIN] << ", " << panels[j][X][MAX] << ", " 
        					<< panels[j][Y][MIN] << ", " << panels[j][Y][MAX] << ", " 	
        					<< panels[j][Z][MIN] << ", " << panels[j][Z][MAX] << ". Dir: "
        					<< dirs[j] << ", " << basisDirs[j] << endl;
        }
        #endif


        //* Combine rows or columns in place if consecutive panels
        //  belong to the same basis function
        if ( (i!=j) && (ind[i]==ind[j]) ){
            ptrP[ ind[i] + nCoefs*(ind[j]-startC[rank]) ] += result*2;
        }else{
            ptrP[ ind[i] + nCoefs*(ind[j]-startC[rank]) ] += result;
        }
    }

    //* Combine sub-matrices to rank0 node
    if( rank==0 ){
        MPI_Status status;
        for ( int i=1; i<numproc; i++ ){
            int copylen = (lastC[i]-startC[i]+1)*nCoefs;
            MPI_Recv(tempP, copylen, MPI_FLOAT, i, 0, MPI_COMM_WORLD, &status);
            double alpha = 1.0f;
            int inc = 1;
            ptrP = dP+(startC[i]*nCoefs);
            daxpy_(&copylen, &alpha, tempP, &inc, ptrP, &inc);
        }
    }else{
        int copylen = (lastC[rank]-startC[rank]+1)*nCoefs;
        MPI_Send(tempP, copylen, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
    }


    delete[] startK;
    delete[] lastK;
    delete[] startC;
    delete[] lastC;
    delete[] tempP;
    delete[] ind;

}



float Caplet::calGalerkinPEntry(int panel1, int panel2){

    //* A nDim-element array of pointers pointing to a nBit-element array
    float *coord_ptr_1[3][4];
    float *coord_ptr_2[3][4];

    if ( this->basisTypes[panel1] == 'F' ){
        int temp = panel1;
        panel1 = panel2;
        panel2 = temp;
    }

    //* Select shapes
    shape_t shape1 = selectShape(panel1);
    shape_t shape2 = selectShape(panel2);

    //* Rotate, mirror, and call proper integrals
    switch( dirs[panel1] ){
    case X:
        switch( basisDirs[panel1] ){
        case Y: // Xy
            this->rotateX2Z(coord_ptr_1, panel1);
            this->rotateX2Z(coord_ptr_2, panel2);
            switch( dirs[panel2] ){
            case X: // Xy_X
                switch( basisDirs[panel2] ){
                case Y: 	// Xy_Xy
                    return intZXZX(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Z: 	// Xy_Xz
                    return intZXZY(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT: 	// Xy_Xf
                    return intZXZF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            case Y: // Xy_Y
                switch( basisDirs[panel2] ){
                case X:		// Xy_Yx
                    return intZXXZ(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Z:		// Xy_Yz
                    return intZXXY(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:	// Xy_Yf
                    return intZXXF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            case Z:
                switch( basisDirs[panel2] ){
                case X:		// Xy_Zx
                    return intZXYZ(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Y: 	// Xy_Zy
                    return intZXYX(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:  // Xy_Zf
                    return intZXYF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            }
        case Z: // Xz
            this->mirrorX2Z(coord_ptr_1, panel1);
            this->mirrorX2Z(coord_ptr_2, panel2);
            switch( dirs[panel2] ){
            case X:
                switch( basisDirs[panel2] ){
                case Y:		// Xz_Xy
                    return intZXZY(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Z:		// Xz_Xz
                    return intZXZX(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:	// Xz_Xf
                    return intZXZF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            case Y:
                switch( basisDirs[panel2] ){
                case X:		// Xz_Yx
                    return intZXYZ(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Z:		// Xz_Yz
                    return intZXYX(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:  // Xz_Yf
                    return intZXYF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
                break;
            case Z:
                switch( basisDirs[panel2] ){
                case X:		// Xz_Zx
                    return intZXXZ(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Y:		// Xz_Zy
                    return intZXXY(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:  // Xz_Zf
                    return intZXXF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            }
        case FLAT: // Xf-?f
            switch( dirs[panel2] ){
            case X: // Xf_Xf
                this->mirrorX2Z(coord_ptr_1, panel1);
                this->mirrorX2Z(coord_ptr_2, panel2);
                return intZFZF(coord_ptr_1, coord_ptr_2);
            case Y: // Xf_Yf
                this->rotateX2Z(coord_ptr_1, panel1);
                this->rotateX2Z(coord_ptr_2, panel2);
                return intZFXF(coord_ptr_1, coord_ptr_2);
            case Z: // Xf_Zf
                this->mirrorX2Z(coord_ptr_1, panel1);
                this->mirrorX2Z(coord_ptr_2, panel2);
                return intZFXF(coord_ptr_1, coord_ptr_2);
            }
        }
    case Y:
        switch( basisDirs[panel1] ){
        case X: 	// Yx
            this->mirrorY2Z(coord_ptr_1, panel1);
            this->mirrorY2Z(coord_ptr_2, panel2);
            switch( dirs[panel2] ){
            case X: // Yx_X
                switch( basisDirs[panel2] ){
                case Y: 	// Yx_Xy
                    return intZXXZ(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Z:		// Yx_Xz
                    return intZXXY(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:  // Yx_Xf
                    return intZXXF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            case Y: // Yx_Y
                switch( basisDirs[panel2] ){
                case X:		// Yx_Yx
                    return intZXZX(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Z:     // Yx_Yz
                    return intZXZY(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:  // Yx_Yf
                    return intZXZF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            case Z: // Yx_Z
                switch( basisDirs[panel2] ){
                case X:		// Yx_Zx
                    return intZXYX(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Y:		// Yx_Zy
                    return intZXYZ(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:  // Yx_Zf
                    return intZXYF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            }
        case Z: 	// Yz
            this->rotateY2Z(coord_ptr_1, panel1);
            this->rotateY2Z(coord_ptr_2, panel2);
            switch( dirs[panel2] ){
            case X: // Yz_X
                switch( basisDirs[panel2] ){
                case Y:		// Yz_Xy
                    return intZXYZ(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Z:		// Yz_Xz
                    return intZXYX(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:	// Yz_Xf
                    return intZXYF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            case Y: // Yz_Y
                switch( basisDirs[panel2] ){
                case X:		// Yz_Yx
                    return intZXZY(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Z:		// Yz_Yz
                    return intZXZX(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:  // Yz_Yf
                    return intZXZF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            case Z: // Yz_Z
                switch( basisDirs[panel2] ){
                case X:		// Yz_Zx
                    return intZXXY(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Y:		// Yz_Zy
                    return intZXXZ(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:  // Yz_Zf
                    return intZXXF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            }
        case FLAT: 	// Yf_?f
            switch( dirs[panel2] ){
            case X:
                this->mirrorY2Z(coord_ptr_1, panel1);
                this->mirrorY2Z(coord_ptr_2, panel2);
                return intZFXF(coord_ptr_1, coord_ptr_2);
            case Y:
                this->mirrorY2Z(coord_ptr_1, panel1);
                this->mirrorY2Z(coord_ptr_2, panel2);
                return intZFZF(coord_ptr_1, coord_ptr_2);
            case Z:
                this->rotateY2Z(coord_ptr_1, panel1);
                this->rotateY2Z(coord_ptr_2, panel2);
                return intZFXF(coord_ptr_1, coord_ptr_2);
            }
        }
    case Z: // Z
        switch( basisDirs[panel1] ){
        case X:		// Zx
            this->rotateZ2Z(coord_ptr_1, panel1);
            this->rotateZ2Z(coord_ptr_2, panel2);
            switch( dirs[panel2] ){
            case X: // Zx_X
                switch( basisDirs[panel2] ){
                case Y:		// Zx_Xy
                    return intZXXY(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Z:		// Zx_Xz
                    return intZXXZ(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:	// Zx_Xf
                    return intZXXF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            case Y:	// Zx_Y
                switch( basisDirs[panel2] ){
                case X:		// Zx_Yx
                    return intZXYX(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Z:		// Zx_Yz
                    return intZXYZ(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:
                    return intZXYF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            case Z: // Zx_Z
                switch( basisDirs[panel2] ){
                case X:		// Zx_Zx
                    return intZXZX(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Y:		// Zx_Zy
                    return intZXZY(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:
                    return intZXZF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            }
        case Y:		// Zy
            this->mirrorY2X(coord_ptr_1, panel1);
            this->mirrorY2X(coord_ptr_2, panel2);
            switch( dirs[panel2] ){
            case X: // Zy_X
                switch( basisDirs[panel2] ){
                case Y:		// Zy_Xy
                    return intZXYX(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Z:		// Zy_Xz
                    return intZXYZ(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:
                    return intZXYF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            case Y: // Zy_Y
                switch( basisDirs[panel2] ){
                case X:		// Zy_Yx
                    return intZXXY(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Z:		// Zy_Yz
                    return intZXXZ(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:
                    return intZXXF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            case Z:	// Zy_Z
                switch( basisDirs[panel2] ){
                case X:		// Zy_Zx
                    return intZXZY(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Y:		// Zy_Zy
                    return intZXZX(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:
                    return intZXZF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            }
        case FLAT:	// Zf_?f
            switch( dirs[panel2] ){
            case X:
                this->rotateZ2Z(coord_ptr_1, panel1);
                this->rotateZ2Z(coord_ptr_2, panel2);
                return intZFXF(coord_ptr_1, coord_ptr_2);
            case Y:
                this->mirrorY2X(coord_ptr_1, panel1);
                this->mirrorY2X(coord_ptr_2, panel2);
                return intZFXF(coord_ptr_1, coord_ptr_2);
            case Z:
                this->rotateZ2Z(coord_ptr_1, panel1);
                this->rotateZ2Z(coord_ptr_2, panel2);
                return intZFZF(coord_ptr_1, coord_ptr_2);
            }
        }
    }

    return 0.0; // dummy return
}


void Caplet::generateRHS(){
	int inc = 1;
	int nRHS = this->nWires*this->nCoefs;
	float alpha = 0.0f;
	sscal_(&nRHS, &alpha, this->coefs, &inc);

	int rowIndex = -1;
	int columnIndex = 0;
	int coefCounter = 0;
	for ( int i=0; i < this->nPanels; i++ ){
		rowIndex 	+= this->indexIncrements[i];
		coefCounter += this->indexIncrements[i];

        //* Dielectric panels: advance rowIndex but contribute nothing to RHS columns
        if ( !this->isDielectric.empty() && this->isDielectric[i] ){
            continue;
        }

		if ( coefCounter > this->nWireCoefs[columnIndex] ){
			coefCounter = 1;
			columnIndex++;
		}

		float (*func)(float x, float w) = &arch;
		switch (this->basisTypes[i]){
		case 'S':
			func = &side;
		case 'A':
			float b = this->panels[i][this->basisDirs[i]][LENGTH];
			float xm = (b)/2;
			float xr = (b)/2;
			int dir = (this->basisDirs[i]+1)%3;
			if ( dir==dirs[i] ){
				dir = (dir+1)%3;
			}
			float w  = this->panels[i][dir][LENGTH];

            //* Gauss quad for computing area integral of shapes
            float val = 0;
            int init_i = 0;
            if ( (gauss_n%2)==1 ){ // odd n
                init_i = 1;
                val += (*gauss::w[gauss_n])[0] * func(xm, w);
            }
            int gauss_n2 = (gauss_n+1)/2;
            for ( int gi = init_i; gi < gauss_n2; gi++ ) {
                float dx = xr * (*gauss::p[gauss_n])[gi];
                val += (*gauss::w[gauss_n])[gi] * ( func(xm+dx, w) + func(xm-dx, w) );
            }
            //* End of Gauss quad

            this->areas[i] = val*xr*w;
		}

		this->coefs[ columnIndex*this->nCoefs + rowIndex ] += this->areas[i];
	}
	scopy_(&nRHS, this->coefs, &inc, this->rhs, &inc);
}


void Caplet::modifyPanelAspectRatio(){
    //* Aspect ratio of all panels cannot exceed MAX_ASPECT_RAIO
	if ( this->isLoaded == false ){
		std::cerr << "ERROR: not yet load any structure " << std::endl;
		std::cerr << "       Unable to check structure validity" << std::endl;

		std::exit(1);
	}
    const float maxRatio = MAX_ASPECT_RATIO;

	std::list<float> 	panelsTemp[3][2];
	std::list<int>		dirsTemp;
	std::list<float>	areasTemp;
	std::list<int>		indexIncrementsTemp;
	std::list<char>		basisTypesTemp;
	std::list<int>		basisDirsTemp;
	std::list<float>	basisZsTemp;
	std::list<float>	basisShiftsTemp;

	int currentWireIndex = 0;
	int currentWirePanelCount = 0;
	int nCurrentWirePanels = this->nWirePanels[0];
	int nFinalPanels = this->nPanels;

	for ( int i=0; i<this->nPanels; i++ ){

		currentWirePanelCount ++;
		if ( currentWirePanelCount > this->nWirePanels[currentWireIndex] ){
			currentWirePanelCount = 1;
			this->nWirePanels[currentWireIndex] = nCurrentWirePanels;
			currentWireIndex ++;
			nCurrentWirePanels = this->nWirePanels[currentWireIndex];
		}
		float ba = this->panels[i][ (this->dirs[i]+1)%3 ][LENGTH] / this->panels[i][ (this->dirs[i]+2)%3 ][LENGTH];
        if (  ( ba > maxRatio )  ||  ( 1/ba > maxRatio )  ){
            //* need to split to reduce aspact ratio of a panel
			int nSplit = 2;
			while ( ( ba/nSplit > maxRatio )  ||  ( 1/ba/nSplit > maxRatio ) ){
				nSplit++;
			}
			int   splitDir = ( ( ba>1 ) ? ( (dirs[i]+1)%3 ) : ( (dirs[i]+2)%3 ) );
			float currentCoord = this->panels[i][splitDir][0];
			float sublength = this->panels[i][splitDir][LENGTH] / nSplit;

			for ( int j=0; j<nSplit; j++){
				panelsTemp[splitDir][0].push_back( currentCoord );
				currentCoord += sublength;
				panelsTemp[splitDir][1].push_back( currentCoord );

				panelsTemp[ (splitDir+1)%3 ][0].push_back( this->panels[i][ (splitDir+1)%3 ][0] );
				panelsTemp[ (splitDir+1)%3 ][1].push_back( this->panels[i][ (splitDir+1)%3 ][1] );
				panelsTemp[ (splitDir+2)%3 ][0].push_back( this->panels[i][ (splitDir+2)%3 ][0] );
				panelsTemp[ (splitDir+2)%3 ][1].push_back( this->panels[i][ (splitDir+2)%3 ][1] );

				dirsTemp.push_back(this->dirs[i]);
				areasTemp.push_back(this->areas[i]/nSplit);
				if ( j==0 ){
                    //* The first sub-panel after split. Needs index increment.
					indexIncrementsTemp.push_back(this->indexIncrements[i]);
				}else{
                    //* The rest of sub-panel after split. No index increment needed.
					indexIncrementsTemp.push_back(0);
				}
				basisTypesTemp.push_back(this->basisTypes[i]);
				basisDirsTemp.push_back(this->basisDirs[i]);
				basisZsTemp.push_back(this->basisZs[i]);

				if ( this->basisDirs[i] == splitDir ){
					if ( this->basisZs[i] > 0 ){ // decaying in the positive direction
						basisShiftsTemp.push_back(this->basisShifts[i] + j*sublength);
					}else{ // decaying in the negative direction
						basisShiftsTemp.push_back(this->basisShifts[i] + (nSplit-1-j)*sublength);
					}
				}else{
					basisShiftsTemp.push_back(this->basisShifts[i]);
				}

				if ( j!=0 ){
					nFinalPanels ++;
					nCurrentWirePanels ++;
				}
			}
        }
        else{
            //* no need to split. simply push back
			panelsTemp[X][0].push_back(this->panels[i][X][0]);
			panelsTemp[X][1].push_back(this->panels[i][X][1]);
			panelsTemp[Y][0].push_back(this->panels[i][Y][0]);
			panelsTemp[Y][1].push_back(this->panels[i][Y][1]);
			panelsTemp[Z][0].push_back(this->panels[i][Z][0]);
			panelsTemp[Z][1].push_back(this->panels[i][Z][1]);

			dirsTemp.push_back(this->dirs[i]);
			areasTemp.push_back(this->areas[i]);
			indexIncrementsTemp.push_back(this->indexIncrements[i]);
			basisTypesTemp.push_back(this->basisTypes[i]);
			basisDirsTemp.push_back(this->basisDirs[i]);
			basisZsTemp.push_back(this->basisZs[i]);
			basisShiftsTemp.push_back(this->basisShifts[i]);
		}
	}
	this->nPanels = nFinalPanels;
	this->nWirePanels[currentWireIndex] = nCurrentWirePanels;

    //* Rebuild panel information based on new panel list with split panels
	delete[] this->panels;
	delete[] this->dirs;
    delete[] this->indexIncrements;
	delete[] this->basisTypes;
	delete[] this->basisDirs;
	delete[] this->basisZs;
	delete[] this->basisShifts;

	this->panels			= new float	[this->nPanels][3][4];
	this->dirs				= new int	[this->nPanels];
	this->areas	 			= new float	[this->nPanels];
	this->indexIncrements 	= new int	[this->nPanels];
	this->basisTypes		= new char	[this->nPanels];
	this->basisDirs			= new int	[this->nPanels];
	this->basisZs			= new float	[this->nPanels];
	this->basisShifts		= new float	[this->nPanels];

	for ( int i=0; i<this->nPanels; i++ ){
		this->panels[i][X][0] = panelsTemp[X][0].front(); panelsTemp[X][0].pop_front();
		this->panels[i][X][1] = panelsTemp[X][1].front(); panelsTemp[X][1].pop_front();
		this->panels[i][Y][0] = panelsTemp[Y][0].front(); panelsTemp[Y][0].pop_front();
		this->panels[i][Y][1] = panelsTemp[Y][1].front(); panelsTemp[Y][1].pop_front();
		this->panels[i][Z][0] = panelsTemp[Z][0].front(); panelsTemp[Z][0].pop_front();
		this->panels[i][Z][1] = panelsTemp[Z][1].front(); panelsTemp[Z][1].pop_front();

		this->panels[i][X][LENGTH] = this->panels[i][X][1] - this->panels[i][X][0];
		this->panels[i][Y][LENGTH] = this->panels[i][Y][1] - this->panels[i][Y][0];
		this->panels[i][Z][LENGTH] = this->panels[i][Z][1] - this->panels[i][Z][0];

		this->panels[i][X][CENTER] = (this->panels[i][X][1] + this->panels[i][X][0])/2;
		this->panels[i][Y][CENTER] = (this->panels[i][Y][1] + this->panels[i][Y][0])/2;
		this->panels[i][Z][CENTER] = (this->panels[i][Z][1] + this->panels[i][Z][0])/2;

		this->dirs[i] = dirsTemp.front(); dirsTemp.pop_front();
		this->areas[i] = areasTemp.front(); areasTemp.pop_front();

		this->indexIncrements[i] = indexIncrementsTemp.front(); indexIncrementsTemp.pop_front();

		this->basisTypes[i] = basisTypesTemp.front(); basisTypesTemp.pop_front();
		this->basisDirs[i] = basisDirsTemp.front(); basisDirsTemp.pop_front();
		this->basisZs[i] = basisZsTemp.front(); basisZsTemp.pop_front();
		this->basisShifts[i] = basisShiftsTemp.front(); basisShiftsTemp.pop_front();
	}
}


shape_t Caplet::selectShape(int panel){
	switch ( this->basisTypes[panel] ){
	case 'A':
		return &arch;
		break;
	case 'S':
		return &side;
		break;
	default:
		return 0;
	}
}


int  Caplet::getNPanels() const{
	return this->nPanels;
}


int  Caplet::getNCoefs() const{
	return this->nCoefs;
}


int  Caplet::getSizeP() const{
	return this->nCoefs * this->nCoefs;
}


int  Caplet::getSizeCoefs() const{
	return this->nCoefs * this->nWires;
}


int  Caplet::getSizeCmat() const{
	return this->nWires * this->nWires;
}


const float* const Caplet::getCmat() const{
	return this->Cmat;
}


float Caplet::compareCmatError(const float* const cmatRef, ERROR_REF option) const{
	using std::abs;

	float* const cmatError = new float[this->getSizeCmat()];
	int n = this->nWires;
	float maxError = 0.0f;
	int maxErrorI = 0;
	int maxErrorJ = 0;

	for ( int i=0; i<n; i++){
		for ( int j=0; j<n; j++ ){

			float errRef = 0.0f;
			switch (option){
			case DIAGONAL:
				errRef = (cmatRef[ i + n*i ]);
				break;
			case SELF:
				errRef = (cmatRef[ i + n*j ]);
			}
			float thisError = abs( (this->Cmat[ i + n*j] - cmatRef[ i + n*j ])/errRef * 100 );
			cmatError[ i + n*j ] = thisError;

			if ( maxError < thisError ){
				maxError  = thisError;
				maxErrorI = i;
				maxErrorJ = j;
			}
		}
	}
	print_matrix(cmatError, n, n, "Cmat Error (%)", 'a');

	std::cout << "Max error = " << std::fixed << maxError << "% at (" << maxErrorI << "," << maxErrorJ << ")" << std::endl;

	switch( option ){
	case DIAGONAL:
		std::cout << "    (w.r.t. the row diagonal)" << std::endl;
		break;
	case SELF:
		std::cout << "    (w.r.t. itself)" << std::endl;
	}
	std::cout << std::endl;

	delete [] cmatError;
	return maxError;
}


float Caplet::compareCmatError(const Caplet *const caplet, ERROR_REF option) const{
	const float* const cmatRef = caplet->getCmat();
	return this->compareCmatError(cmatRef, option);
}


float Caplet::compareCmatError(const std::string filename, ERROR_REF option) const{
	const int n = this->nWires;

	float* cmatRef = new float[n*n];

	std::ifstream ifile(filename.c_str());
	if (!ifile){
		std::cerr << "ERROR: cannot open the file: " << filename << std::endl;
		std::exit(1);
	}

	std::string		lineTemp;
	for ( int i=0; i<n; i++ ){
		getline(ifile, lineTemp);
		std::stringstream	stringTokenizer(lineTemp);
		for (int j=0; j<n; j++){
			stringTokenizer >> cmatRef[ i + n*j ];
		}
	}

	float maxError = this->compareCmatError(cmatRef, option);
	delete [] cmatRef;

	return maxError;
}



//__________________________________________________________
//*
//* PRINT UTILITIES
//*
void Caplet::printP(std::ostream &out){
    if ( this->isLoaded && this->isSolved ){
        switch(this->mode){
        case FAST_GALERKIN:
            print_matrix(this->P, this->nCoefs, this->nCoefs, "", 'u', out);
            break;
        case DOUBLE_GALERKIN:
        case DOUBLE_COLLOCATION:
            print_matrix(this->dP, this->nCoefs, this->nCoefs, "", 'u', out);
            cout << "here" << endl;
            cout << dP[0] << ", " << dP[1] << endl;
            out << "test" << endl;
            break;
        default:
            ;
        }
    }
}


void Caplet::printCoefs(){
    if ( this->isLoaded && this->isSolved ){
        switch(this->mode){
        case FAST_GALERKIN:
            print_matrix(this->coefs, this->nCoefs, this->nWires, "");
            break;
        case DOUBLE_GALERKIN:
        case DOUBLE_COLLOCATION:
            print_matrix(this->dcoefs, this->nCoefs, this->nWires, "");
            break;
        default:
            ;
        }
    }
}


void Caplet::printRHS(std::ostream &out){
    if ( this->isLoaded && this->isSolved ){
        switch(this->mode){
        case FAST_GALERKIN:
            print_matrix(this->rhs, this->nCoefs, this->nWires, "", 'a', out);
            break;
        case DOUBLE_GALERKIN:
        case DOUBLE_COLLOCATION:
            print_matrix(this->drhs, this->nCoefs, this->nWires, "", 'a', out);
            break;
        default:
            ;
        }
    }
}


void Caplet::printCmat(){
    if ( this->isLoaded && this->isSolved ){
        switch(this->mode){
        case FAST_GALERKIN:
            print_matrix(this->Cmat, this->nWires, this->nWires, "");
            break;
        case DOUBLE_GALERKIN:
        case DOUBLE_COLLOCATION:
            print_matrix(this->dCmat, this->nWires, this->nWires, "");
            break;
        default:
            ;
        }
    }
}


void Caplet::printPanel(int index){
    if (isLoaded == false){
        return;
    }
    if (index >= nPanels){
        return;
    }
    for (int i=0; i<3; ++i){
        for (int j=0; j<4; ++j){
            cout << setw(12) << panels[index][i][j];
        }
        cout << endl;
    }
}



//__________________________________________________________
//*
//* DEBUG UTILITIES
//*
bool Caplet::isPanelAspectRatioValid(){
    if ( this->isLoaded == false ){
        std::cerr << "ERROR: not yet load any structure " << std::endl;
        std::cerr << "       Unable to check structure validity" << std::endl;

        return false;
    }
    const float maxRatio = 51;

    for ( int i=0; i< this->nPanels; i++ ){
        float ratio = panels[i][(dirs[i]+1)%3][LENGTH] / panels[i][(dirs[i]+2)%3][LENGTH];

        if ( ratio > maxRatio || 1/ratio > maxRatio ){
            std::cerr << "ERROR: panel ratio check fails" << std::endl;
            std::cerr << "    Panel: " << i << std::endl;
            std::cerr << "    x:" << panels[i][X][0] << ", " << panels[i][X][1] << ", " << panels[i][X][2] << ", " << panels[i][X][3] << std::endl;
            std::cerr << "    y:" << panels[i][Y][0] << ", " << panels[i][Y][1] << ", " << panels[i][Y][2] << ", " << panels[i][Y][3] << std::endl;
            std::cerr << "    z:" << panels[i][Z][0] << ", " << panels[i][Z][1] << ", " << panels[i][Z][2] << ", " << panels[i][Z][3] << std::endl;

            return false;
        }
    }
    return true;
}



} // end of namespace caplet