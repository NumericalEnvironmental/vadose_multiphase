// MFP - a finite-difference code for modeling multiphase flow in 1-D, 2-D, or 3-D porous media
// model solves PDE for total fluid pressure, then calculates fluxes of separate phases over given time step before updating

import std.stdio; 			// I/O and file system
import std.math; 			// math function library
import std.algorithm; 		// algorithm for working with ranges, sequences, etc. 'Sorting' is a submodule, which is used herein.
import std.array; 			// array operations
import std.string; 			// string function support
import std.conv; 			// automatic conversions between types
import std.typecons; 		// misc. meta-programming support, including type constructors, tuples, etc.


// classes ...

class Liquid{

	// properties of liquids (i.e., water and NAPL)
	string name;
	double rho,u,beta_nw,beta_na;
	
	this(string name,double rho,double u,double beta_nw,double beta_na){
		this.name = name;			// name
		this.rho = rho; 			// density
		this.u = u; 				// viscosity
		this.beta_nw = beta_nw;		// 1 - air-water/NAPL-water surface tension ratio
		this.beta_na = beta_na;		// air-water/NAPL-air surface tension ratio
		}
	
	}

class Connection{

	// cell-to-cell connection properties
	int icell_0,icell_1,dir;
	double A,d;

	this(int icell_0,int icell_1,int dir,double A,double d){
		this.icell_0 = icell_0;
		this.icell_1 = icell_1;
		this.dir = dir;
		this.A = A;
		this.d = d;
		}
	
	double C(int phase,Cell[] cell,Liquid[] liquid,double g){
		// compute fluid conductance between cell_1 and cell_0, using upstream weighting
		// if grad < 0, flow is from cell 1 into cell 0
		double K_upstream;
		if (Grad(phase,cell,liquid,g) < 0.0){K_upstream = cell[icell_1].kr(phase)*cell[icell_1].k;}
		else{K_upstream = cell[icell_0].kr(phase)*cell[icell_0].k;}
		return (A/d)*K_upstream/liquid[phase].u;
		}
	
	double Grad(int phase,Cell[] cell,Liquid[] liquid,double g){
		// potential gradient from cell_0 into cell_1
		return (cell[icell_0].P-cell[icell_1].P) + liquid[phase].rho*g*(cell[icell_0].z-cell[icell_1].z);
		}	
	
	double Flux(int phase,Cell[] cell,Liquid[] liquid,double g){
		// phase-specific fluid flux from cell_0 into cell_1
		return C(phase,cell,liquid,g) * Grad(phase,cell,liquid,g);
		}
	
	}
	
class Source{

	// source terms (per well, per specified period)
	int i;
	int[] icell;
	double[] Q_t;
	double[] t0;
	double[] tf;
	double[] Q;
	string line_input;
	string[] parsed_input;	
	
	this(int num_cells,string file_name){
		this.Q = ZeroesVector(num_cells); 		// assign background flux per cell, as a placeholder
		// read input file
		auto input_file = File(file_name,"r");
		while (!input_file.eof()) {
			line_input = input_file.readln();
			if (i > 0) { 				// ignore header line in input file 
				parsed_input = split(line_input);
				this.icell ~= to!int(parsed_input[0]);
				this.Q_t ~= to!double(parsed_input[1]); 		// flux into icell, starting at time t0 and ending at time tf
				this.t0 ~= to!double(parsed_input[2]);				
				this.tf ~= to!double(parsed_input[3]);				
				}
			++i;
			}
		input_file.close();		
		}
		
	double[] Q_time(double time){
		// return array of Q values, which will include those for cells with finite source terms at this time
		for (int i = 0; i < icell.length; ++i){
			if (time >= t0[i] && time <= tf[i]){Q[icell[i]] = Q_t[i];}
			}
		return Q;
		}
	
	}	
	
class Grid{

	// Grid: sizing of grid, interfacial areas, etc.
	double x0,xf,y0,yf,z0,zf,Ax,Ay,Az,dx,dy,dz,vol;
	int i,nx,ny,nz,N;
	string line_input;
	string[] parsed_input;
	
	this(){
	
		// read grid settings
		auto input_file = File("grid.txt","r");
		while (!input_file.eof()) {
			line_input = input_file.readln();
			if (i > 0) { 				// ignore header line in input file 
				parsed_input = split(line_input);
				switch (parsed_input[0]) {
					case "start": 					
						this.x0 = to!double(parsed_input[1]);
						this.y0 = to!double(parsed_input[2]);
						this.z0 = to!double(parsed_input[3]);						
						break;
					case "end": 		 			
						this.xf = to!double(parsed_input[1]);
						this.yf = to!double(parsed_input[2]);
						this.zf = to!double(parsed_input[3]);	
						break;
					default: 						
						assert(parsed_input[0] == "N");
						this.nx = to!int(parsed_input[1]);
						this.ny = to!int(parsed_input[2]);
						this.nz = to!int(parsed_input[3]);
						break;
					}						
				}
			++i;
			}		
		input_file.close();
		
		// assign default geometry
		this.N = nx*ny*nz;
		this.dx = (xf-x0)/nx;
		this.dy = (yf-y0)/ny;
		this.dz = (zf-z0)/nz;
		this.Ax = dy*dz;
		this.Ay = dx*dz;
		this.Az = dx*dy;
		this.vol = dx*dy*dz;
		}

	}

class Model{
	
	// Model: (1) intercell fluxes, (2) matrix assembly from cell object collection, and (3) solution of linear equation system using SOR
	static double gamma = 0.5; 		// finite-difference time weighting factor (0.5 < gamma <= 1.0)
	double g,dt_0,dt_min,dt_max,dP_max,ssv,dt_down_f,dt_up_f,t_end,print_interval;
	int num_cells,solver_type;
	string line_input;	
	string[] parsed_input;
	
	this(int num_cells,double g) {
		this.num_cells = num_cells;
		this.g = g;							// gravitational acceleration (included in class)
		this.gamma = gamma;
		// read model controls file
		auto input_file = File("controls.txt","r");
		while (!input_file.eof()) {
			line_input = input_file.readln();
			parsed_input = split(line_input);
			switch (parsed_input[0]) {
				case "dt_0": 						// starting time step size
					this.dt_0 = to!double(parsed_input[1]);
					break;
				case "dt_min": 		 				// minimum time step size
					this.dt_min = to!double(parsed_input[1]);
					break;
				case "dt_max": 						// maximum time step
					this.dt_max = to!double(parsed_input[1]);
					break;
				case "dP_max": 		 				// maximum implied pressure change per time step
					this.dP_max = to!double(parsed_input[1]);
					break;		
				case "steady-state_dSdt": 		 				// steady-state threshold criteria (dS/dt)
					this.ssv = to!double(parsed_input[1]);
					break;
				case "dt_down_f": 					// factor used to reduce time step size if dP_max is exceeded
					this.dt_down_f = to!double(parsed_input[1]);
					break;
				case "dt_up_f": 		 			// factor used to grow time step, when dP_max is not exceeded, until dt_max is reached
					this.dt_up_f = to!double(parsed_input[1]);
					break;					
				case "t_end": 		 				// purely saturated condition; compute specific storage pressure response
					this.t_end = to!double(parsed_input[1]);
					break;		
				case "print_interval": 		 				// print output time interval
					this.print_interval = to!double(parsed_input[1]);
					break;
				default: 							// matrix solver type (0 = successive over-relaxation; 1 = tridiagonal)
					assert(parsed_input[0] == "solver_type");
					this.solver_type = to!int(parsed_input[1]);
					break;
				}		
			}		
		input_file.close();
		}

	void Simulate(Cell[] cell,Connection[] connection,Liquid[] liquid,Source source_w,Source source_n,double time){
		// step through simulation ... 
		double dt,print_time,dt_preprint,dSw_min,dSw_max,dSn_min,dSn_max;
		double[] dP; 										// changes in fluid pressure implied by matrix solver
		double[] dSw;
		double[] dSn;
		double[] dP_s;
		double[] dP_s0;
		double[] Qw;
		double[] Qn;
		string file_name;
		bool i_pass,print_flag;
		const double del_t = 1e-7; 							// tolerance for difference between 'time' and 'print time'
		dSw = ZeroesVector(num_cells);
		dSn = ZeroesVector(num_cells);		
		dP_s = ZeroesVector(num_cells);
		dP_s0 = ZeroesVector(num_cells);			
		// write initial conditions & fluxes to files
		for (int i = 0; i < num_cells; ++i){cell[i].P = cell[i].Pressure(liquid,g);} 	// set initial fluid pressure as a function of saturation(s)
		file_name = "init_condts.txt";
		WriteState(file_name,cell,liquid,g);
		file_name = "init_velocities.txt";		
		WriteVelocities(file_name,cell,connection,liquid,g);
		WriteFluxes(connection,cell,liquid,g);
		dt = dt_0;		
		do {
			print_time = time + (print_interval - time % print_interval); 		// next print output step; adjust dt if needed
			dt_preprint = dt;			
			if (time + dt > print_time){
				dt = print_time - time;
				print_flag = 1;
				}			
			else{print_flag = 0;}
			do {							// iterate until all dP ,= dP_max or dt = dt_min; default assumption is that this condition is true (i_pass == 1, until corrected)
				i_pass = 1;
				Qw = source_w.Q_time(time); 				// return source fluxes at beginning of current time step
				Qn = source_n.Q_time(time);		
				dP = SolveStep(cell,connection,liquid,Qw,Qn,dt,dP_s,dP_s0);		// provisional total saturation change
				for (int i = 0; i < num_cells; ++i){		// check implied pressure change
					if (abs(dP[i]) > dP_max){
						// time step is too large for at least one cell
						i_pass = 0;
						dt *= dt_down_f;
						assert(dt>=dt_min,"dt too small at cell " ~to!string(i) ~" where P = " ~to!string(cell[i].P) ~" and dP = " ~to!string(dP[i]));
						break;
						}
					}
				} while (i_pass == 0);
			// implied changes in saturation
			dSw = DeltaSat(0,dP,cell,connection,liquid,num_cells,Qw,g,dt,gamma);
			dSn = DeltaSat(1,dP,cell,connection,liquid,num_cells,Qn,g,dt,gamma);
			//update cell states
			dSw_min = dSw[0];
			dSw_max = dSw[0];
			dSn_min = dSn[0];
			dSn_max = dSn[0];
			for (int i = 0; i < num_cells; ++i) { 	
				cell[i].Sw += dSw[i];
				cell[i].Sn += dSn[i];			
				cell[i].P = cell[i].Pressure(liquid,g);
				dSw_min = min(dSw_min,dSw[i]);
				dSw_max = max(dSw_max,dSw[i]);
				dSn_min = min(dSn_min,dSn[i]);
				dSn_max = max(dSn_max,dSn[i]);				
				}
			time += dt;
			// check if steady-state criteria met
			if (abs(dSw_min/dt) <= ssv && abs(dSw_max/dt) <= ssv && abs(dSn_min/dt) <= ssv && abs(dSn_max/dt) <= ssv){
				file_name = "state_" ~ to!string(time) ~ ".txt";
				WriteState(file_name,cell,liquid,g);
				file_name = "velocities_" ~ to!string(time) ~ ".txt";				
				WriteVelocities(file_name,cell,connection,liquid,g);
				writeln("Steady-state reached by time = ",to!string(time),"; dt = ",to!string(dt_preprint));
				writeln("\tdSw/dt_min = ",to!string(dSw_min/dt),"; dSw/dt_max = ",to!string(dSw_max/dt));
				writeln("\tdSn/dt_min = ",to!string(dSn_min/dt),"; dSn/dt_max = ",to!string(dSn_max/dt));				
				return;
				}
			// write output, as warranted
			if (abs(time-print_time) <= del_t){
				file_name = "state_" ~ to!string(time) ~ ".txt";
				WriteState(file_name,cell,liquid,g);
				file_name = "velocities_" ~ to!string(time) ~ ".txt";				
				WriteVelocities(file_name,cell,connection,liquid,g);
				writeln("Wrote output at time = ",to!string(time),"; dt = ",to!string(dt_preprint));
				writeln("\tdSw/dt_min = ",to!string(dSw_min/dt),"; dSw/dt_max = ",to!string(dSw_max/dt));
				writeln("\tdSn/dt_min = ",to!string(dSn_min/dt),"; dSn/dt_max = ",to!string(dSn_max/dt));				
				}
			// adjust time step, as warranted
			if(print_flag == 1){dt = dt_preprint;}
			else{dt = min(dt*dt_up_f,dt_max);}
			} while (time < t_end);
		}
		
	double[] SolveStep(Cell[] cell,Connection[] connection,Liquid[] liquid,double[] Qw,double[] Qn,double dt,double[] dP,double[] dP0){
		// set up and solve matrix of total flow equations (both phases, summed together) for a given time step
		double[][] a;			
		double[] a_d,b;
		int iternum = 0;
		bool iterSOR = 1;	
		auto submatrix = BuildMatrix(cell,connection,liquid,Qw,Qn,dt);
		a_d = submatrix[0];
		a = submatrix[1];
		b = submatrix[2];
		if (solver_type==0){
			while (iterSOR == 1) {
				++iternum;
				dP = SOR(a_d,a,b,dP,dP0,cell); 								// iterate on solution
				iterSOR = CheckConvgSOR(dP,dP0,iternum);                       	// check for convergence
				for (int i = 0; i < num_cells; ++i) {dP0[i] = dP[i];}
				}	
			}
		else{dP = Tridiagonal(a_d,a,b,cell);}
		return dP;
		}
		
	Tuple!(double[],double[][],double[]) BuildMatrix(Cell[] cell,Connection[] connection,Liquid[] liquid,double[] Qw,double[] Qn,double dt) {
		double conduct;
		int k;
		double[][] a;					// off-diagonal elements
		double[] a_d; 					// diagonal elements
		double[] b; 						// explicit (known) vector
		a.length = num_cells;
		for (int i = 0; i < num_cells; ++i) {
			a_d ~= cell[i].dSdP(liquid,g)*cell[i].phi*cell[i].vol/dt; 		
			b ~= Qw[i]+Qn[i];
			}			
		for (int i = 0; i < num_cells; ++i){
			k = 0;
			foreach(j; cell[i].c_index){
				conduct = connection[j].C(0,cell,liquid,g) + connection[j].C(1,cell,liquid,g); 						// total fluid conductivity (both phases, summed)
				a_d[i] += gamma * conduct;
				a[i] ~= -gamma * conduct;
				b[i] -= -cell[i].dir_index[k] * (connection[j].Flux(0,cell,liquid,g) - connection[j].Flux(1,cell,liquid,g));
				++k;
				}
			}
		return tuple(a_d,a,b);
		}
		
	double[] Tridiagonal(double[] b,double[][] a_od,double[] r,Cell[] cell){
		// tri-diagonal banded matrix solver; translated from Numerical Recipes in Fortran, Press et al. (1992); variable nomenclature maintained 
		double bet;
		double[] a;
		double[] c;
		double[] u;
		double[] gam;
		a = ZeroesVector(num_cells);
		c = ZeroesVector(num_cells);
		u = ZeroesVector(num_cells);   
		gam = ZeroesVector(num_cells);    
		for (int i = 0; i < num_cells-1; ++i){
			c[i] = a_od[i][sgn(i)];					// for a 1-D column, a_od[i][1] is the second (downstream) connecting cell in row i (except for 1st cell)
			a[i+1] = a_od[i+1][0];					// a_od[i][1] is the first (upstream) connecting cell in row i+1
			}
		bet = b[0];
		u[0] = r[0]/bet;
		for (int j = 1; j < num_cells-1; ++j){		// decomposition and forward substitution
			gam[j] = c[j-1]/bet;
			bet = b[j] - a[j]*gam[j];
			u[j] = (r[j] - a[j]*u[j-1])/bet;
			}
		for (int j = num_cells-2; j >= 0; --j){ 	// back substitution
			u[j] = u[j] - gam[j+1]*u[j+1];
			}
		return u;
		}	
		
	double[] SOR(double[] a_d,double[][] a,double[] b,double[] x,double[] x0,Cell[] cell) {
		// solve flow equation matrix by successive over-relaxation (SOR) method	
		const double omega = 0.5; 	// SOR relaxation factor
		double sum_0,sum_1;
		for (int i = 0; i < num_cells; ++i) {
			sum_0 = 0.0;
			sum_1 = 0.0;
			foreach(int k, j; cell[i].n_index) {
				if (j < i) {sum_0 += a[i][k] * x[j];}
				else {sum_1 += a[i][k] * x0[j];}
				}
			x[i] = (1.0 - omega) * x0[i] + (omega/a_d[i]) * (b[i] - sum_0 - sum_1);			
			}
		return x;
		}	

	bool CheckConvgSOR(double[] dP,double[] dP0,int iternum) {
		// check convergence of SOR procedure	
		const int maxiterSOR = 500;
		const double relconvgSOR = 1e-12; 		// dP correction per iteration is to be less than or equal to this value everywhere
		const double dP_threshold = 1e-5;
		bool iterSOR = 0;
		if (iternum == 1) { 				
			// force at least two iterations to occur (so that dP_0 array can be populated for comparison)
			iterSOR = 1;
			return iterSOR;
			}
		if (iternum > maxiterSOR) {
			writeln("\t \t \t WARNING: SOR iteration has exceeded max_iter_SOR.");
			iterSOR = 0;
			return iterSOR;
			}
		for (int i = 0; i < num_cells; ++i) {
			if (abs(dP[i]) > dP_threshold) {
				if (abs(dP[i]-dP0[i])/abs(dP[i]) > relconvgSOR) {
					// convergence has not occurred for this cell
					iterSOR = 1;
					break;
					}
				}
			}
		return iterSOR;
		}
	}

class Cell{

	// Cell: specific attributes and methods that act upon individual cells
	double x,y,z,Sw,Sn,k,phi,vol,alpha,n,m,Swr,Snr,Snr_m,Ss,P,g;
	int[] c_index;
	int[] n_index;
	int[] dir_index;
	this(double x,double y,double z,double Sw,double Sn,double k,double phi,double vol,double alpha,double n,double Swr,double Snr,double Snr_m,double Ss,double g){
		this.x = x;
		this.y = y;
		this.z = z;
		this.c_index = []; 				// array of connections linked to this cell
		this.n_index = []; 				// corresponding cell indices associated with each connection
		this.dir_index = [];			// associated connection direction (+ = flow into cell, - = flow out of cell)
		this.Sw = Sw;
		this.Sn = Sn;
		this.k = k;
		this.phi = phi;					// porosity
		this.vol = vol;					// cell volume; normally implied by grid geometry but can be set to a large value to effect fixed pressure and saturation
		this.alpha = alpha; 			// Van Genuchten parameters
		this.n = n;
		this.m = 1.0 - 1.0/n;
		this.Swr = Swr;
		this.Snr = Snr;
		this.Snr_m = Snr_m;
		this.Ss = Ss;
		this.P = 0.0;					// fluid pressure, based on saturation (placeholder value)
		}
	
	double kr(int phase){
		// relative permeability as a function of capillary pressure head
		double rel_perm;
		double Sw_e,Sn_e,St;		
		auto Se = SatEff();
		Sw_e = Se[0];
		Sn_e = Se[1];
		St = Se[2];		
		if (phase == 0){rel_perm = min(Sw_e,1.0)^^0.5 * (1.0 - (1.0-min(Sw_e,1.0)^^(1.0/m))^^m)^^2.0;} 		// water
		else{																								// NAPL
			if (Sn > Snr_m){rel_perm = (min(St,1.0)-min(Sw_e,1.0))^^0.5 * ((1.0-min(Sw_e,1.0)^^(1.0/m))^^m - (1.0-min(St,1.0)  ^^(1.0/m))^^m)^^2.0;}
			else{rel_perm = 0.0;} 																			// immobile NAPL
			}
		return rel_perm;		
		}
	
	double AlphaEff(Liquid[] liquid){
		// calculate a weighted Van Genuchten alpha for a binary liquid mixture, possibly including an air phase as well
		const double epsilon = 1e-6; 		// factor used to prevent divide-by-zero errors for calculations yielding ~infinity
		double Sa,aw,nw,na,aw0,nw0,na0,norm,a;
		Sa = max(1.0 - min(Sw,1.0) - min(Sn,1.0),0.0); 			
		aw = min(Sw,1.0) + Sa;
		nw = min(Sw,1.0) + min(Sn,1.0);
		na = min(Sn,1.0) + Sa;
		aw0 = 1.0/abs(1.0 - aw + epsilon);
		nw0 = 1.0/abs(1.0 - nw + epsilon);
		na0 = 1.0/abs(1.0 - na + epsilon);
		norm = aw0 + nw0 + na0;
		a = aw0/norm*alpha + nw0/norm*liquid[1].beta_nw*alpha + na0/norm*liquid[1].beta_na*alpha;
		return a;
		}
	
	Tuple!(double,double,double) SatEff(){
		// effective liquid and total liquid saturation (note: can be slightly > 1 under saturated conditions)
		double Sw_e,Sn_e,St;
		Sw_e = (Sw - Swr)/(1.0 - Swr);
		Sn_e = (Sn - Snr)/(1.0 - Snr);
		St = (Sw + Sn - Swr - Snr)/(1.0 - Swr - Snr);
		return tuple(Sw_e,Sn_e,St);
		}

	double Pressure(Liquid[] liquid,double g){
		// fluid pressure as a function of saturation
		// delS = putative total liquid saturation change; passed to method to test potential pressure response
		double Sw_e,Sn_e,St,delSt,alpha_current,avg_rho,Pf;		
		auto Se = SatEff();
		Sw_e = Se[0];
		Sn_e = Se[1];
		St = Se[2];		
		alpha_current = AlphaEff(liquid);
		avg_rho = Sw*liquid[0].rho/(Sw + Sn) + Sn*liquid[1].rho/(Sw + Sn);
		if (Sw + Sn < 1.0){Pf = -exp(log(exp(-log(St)/m)-1.0)/n - log(alpha_current));}
		else{Pf = (Sw + Sn - 1.0)*phi*avg_rho*g/Ss;}
		return Pf;
		}
		
	double dSdP(Liquid[] liquid,double g){
		// change in saturation per change in fluid pressure
		double alpha_current,avg_rho,S_prime_num,S_prime;
		alpha_current = AlphaEff(liquid);		
		avg_rho = Sw*liquid[0].rho/(Sw + Sn) + Sn*liquid[1].rho/(Sw + Sn);
		if (Sw + Sn < 1.0){ 																							
			// unsaturated condition
			S_prime_num = -(alpha_current*abs(P))^^n * ((alpha_current*abs(P))^^n + 1.0)^^(-m-1.0) * m * n/(P);
			S_prime = S_prime_num/(1.0 - Swr - Snr);
			}
		else{S_prime = Ss/(phi*avg_rho*g);}			// saturated condition
		return S_prime;
		}

	}
	
// utility functions ...


double Average(double[] a){
	// return the arithmetic mean
	double sum = 0.0;
	foreach (b; a) {sum += b;}
	return sum/a.length;
	}	

double[] ZeroesVector(int M) {
	// return a vector of doubles of length M, all initialized to zero
	double[] a;
	for (int i = 0; i < M; ++i) {a ~= 0.0;}
	return a;
	}

	
// set up functions (only used once, at the beginning of the simulation) ...
	

Cell[] CellIndices(Cell[] cell,Connection[] connection){
	// populate cell connection index lists
	int[][] row;
	row.length = cell.length;
	for (int i = 0; i < connection.length; ++i){
		// raw connection lists (unsorted)
		row[connection[i].icell_0] ~= i;
		row[connection[i].icell_1] ~= i;		
		}
	for (int i = 0; i < cell.length; ++i) {
		cell[i].c_index = row[i].sort; 				// sort connection indices for each cell (to configure for SOR matrix solver)
		// populate connecting cell lists
		foreach(j; cell[i].c_index){
			// j = index number for connection linked to cell i
			if (connection[j].icell_0 == i){
				cell[i].n_index ~= connection[j].icell_1;
				cell[i].dir_index ~= -1;
				}
			else {
				cell[i].n_index ~= connection[j].icell_0;
				cell[i].dir_index ~= 1;				
				}
			}
		}
	return cell;
	}	
	
Connection[] Connect(Grid grid){
	Connection[] connection; 		// declare 'connection' as a dynamic array of Connection objects
	// cell connections parallel to x-axis
	int icell_0,icell_1;
	for (int k = 0; k < grid.nz; ++k) {
		for (int j = 0; j < grid.ny; ++j) {
			for (int i = 0; i < grid.nx-1; ++i) {
				icell_0 = k*grid.nx*grid.ny + j*grid.nx + i; 
				icell_1 = k*grid.nx*grid.ny + j*grid.nx + i + 1;
				Connection connection_member = new Connection(icell_0,icell_1,0,grid.Ax,grid.dx);
				connection ~= connection_member;
				}
			}
		}
	// cell connections parallel to y-axis
	for (int k = 0; k < grid.nz; ++k) {
		for (int j = 0; j < grid.ny-1; ++j) {
			for (int i = 0; i < grid.nx; ++i) {
				icell_0 = k*grid.nx*grid.ny + i + j*grid.nx; 
				icell_1 = k*grid.nx*grid.ny + i + j*grid.nx + grid.nx;
				Connection connection_member = new Connection(icell_0,icell_1,1,grid.Ay,grid.dy);
				connection ~= connection_member;
				}
			}
		}
	// cell connections parallel to z-axis
	for (int k = 0; k < grid.nz-1; ++k) {
		for (int j = 0; j < grid.ny; ++j) {
			for (int i = 0; i < grid.nx; ++i) {
				icell_0 = k*grid.nx*grid.ny + j*grid.nx + i;
				icell_1 = (k+1)*grid.nx*grid.ny + j*grid.nx + i;
				Connection connection_member = new Connection(icell_0,icell_1,2,grid.Az,grid.dz);
				connection ~= connection_member;
				}
			}
		}
	return connection;
	}	
	
Cell[] ConstructCells(Grid grid,Liquid[] liquid,double g){
	// populate cell object set
	int line_index,cell_index;
	double Sw,Sn,perm,phi,vol,alpha,n,Swr,Snr,Snr_m,Ss;
	double x[];
	double y[];
	double z[];
	string line_input;
	string[] parsed_input;
	Cell[] cell; 		// declare 'cell' as a dynamic array of Cell objects
	// assign coordinates
	for (int k = 0; k < grid.nz; ++k) {
		for (int j = 0; j < grid.ny; ++j) {
			for (int i = 0; i < grid.nx; ++i) {
				x ~= i*grid.dx + 0.5*grid.dx; 
				y ~= j*grid.dy + 0.5*grid.dy;
				z ~= k*grid.dz + 0.5*grid.dz;
				}
			}
		}	
	// read cell properties file
	auto input_file = File("cells.txt","r");
	while (!input_file.eof()) {
		line_input = input_file.readln();
		if (line_index > 0) { 							// ignore header line in input file 
			parsed_input = split(line_input);		
			if (parsed_input[0]=="default") {
				Sw = to!double(parsed_input[1]);
				Sn = to!double(parsed_input[2]);			
				perm = to!double(parsed_input[3]);
				phi = to!double(parsed_input[4]);
				vol = grid.vol;
				alpha = to!double(parsed_input[6]); 			// note: volume not read in 
				n = to!double(parsed_input[7]);
				Swr = to!double(parsed_input[8]);
				Snr = to!double(parsed_input[9]);
				Snr_m = to!double(parsed_input[10]);				
				Ss = to!double(parsed_input[11]);
				for (int i = 0; i < grid.N; ++i) {
					Cell cell_member = new Cell(x[i],y[i],z[i],Sw,Sn,perm,phi,vol,alpha,n,Swr,Snr,Snr_m,Ss,g);
					cell ~= cell_member;			
					}
				}
			else {
				// special cell, with properties unlike default
				cell_index = to!int(parsed_input[0]);
				cell[cell_index].Sw = to!double(parsed_input[1]);
				cell[cell_index].Sn = to!double(parsed_input[2]);
				cell[cell_index].k = to!double(parsed_input[3]);
				cell[cell_index].phi = to!double(parsed_input[4]);
				cell[cell_index].vol = to!double(parsed_input[5]);
				cell[cell_index].alpha = to!double(parsed_input[6]);
				cell[cell_index].n = to!double(parsed_input[7]);
				cell[cell_index].Swr = to!double(parsed_input[8]);
				cell[cell_index].Snr = to!double(parsed_input[9]);
				cell[cell_index].Snr_m = to!double(parsed_input[10]);				
				cell[cell_index].Ss = to!double(parsed_input[11]);
				cell[cell_index].m = 1.0 - 1.0/cell[cell_index].n; 		// internal cell property modifications/placeholders
				cell[cell_index].P = 0.0;	
				}
			}
		++line_index;
		}		
	input_file.close();		
	return cell;
	}

Liquid[] LiquidProps(){
	// read in liquid properties and return 
	int line_index;
	double rho,u,beta_nw,beta_na;
	string line_input,name;
	string[] parsed_input;
	Liquid[] liquid;
	auto input_file = File("liquids.txt","r");
	while (!input_file.eof()) {
		line_input = input_file.readln();
		if (line_index > 0) { 				// ignore header line in input file 
			parsed_input = split(line_input);
			name = parsed_input[0];	
			rho = to!double(parsed_input[1]);
			u = to!double(parsed_input[2]);	
			beta_nw = to!double(parsed_input[3]);				
			beta_na = to!double(parsed_input[4]);
			Liquid liquid_member = new Liquid(name,rho,u,beta_nw,beta_na);
			liquid ~= liquid_member;			
			}
		++line_index;
		}
	input_file.close();		
	return liquid;
	}
	
double[] DeltaSat(int phase,double[] dP,Cell[] cell,Connection[] connection,Liquid[] liquid,int num_cells,double[] Q,double g,double dt,double gamma){
	// compute saturation changes across the cell set, per phase
	double[] dSp;
	double[] sum_term;
	double flux;
	int i,j;
	dSp = ZeroesVector(num_cells);
	sum_term = ZeroesVector(num_cells);
	for (int iconnect = 0; iconnect < connection.length; ++iconnect){
		i = connection[iconnect].icell_0;
		j = connection[iconnect].icell_1;
		flux = connection[iconnect].C(phase,cell,liquid,g) * (((cell[j].P + gamma*dP[j]) - (cell[i].P + gamma*dP[i])) + liquid[phase].rho*g*(cell[j].z-cell[i].z));
		sum_term[i] += flux;
		sum_term[j] -= flux;
		}
	for (int k = 0; k < num_cells; ++k){dSp[k] = dt*(Q[k] + sum_term[k])/(cell[k].phi*cell[k].vol);}
	return dSp;
	}
	
void WriteState(string file_name,Cell[] cell,Liquid[] liquid,double g) {
	// write results to output file
	auto output_file = File(file_name, "w");
	output_file.writeln("cell","\t","x","\t","y","\t","z","\t","Sw","\t","Sn","\t","P","\t","dSdP");
	for (int i = 0; i < cell.length; ++i) {output_file.writeln(i,"\t",cell[i].x,"\t",cell[i].y,"\t",cell[i].z,"\t",cell[i].Sw,"\t",cell[i].Sn,"\t",cell[i].P,"\t",cell[i].dSdP(liquid,g));}
	output_file.close();
	}	

Tuple!(double,double,double) Velocity(int cell_index,int phase,Connection[] connection,Cell[] cell,Liquid[] liquid,double g){
	// summarize current velocity vector for 'phase' for cell[cell_index], averaged across adjoining cell interfaces
	double vx = 0.0;
	double vy = 0.0;
	double vz = 0.0;
	int vx_count,vy_count,vz_count;
	foreach(j; cell[cell_index].c_index){
		switch (connection[j].dir) {
			case 0: 						// x-axis connection
				vx += connection[j].Flux(phase,cell,liquid,g)/connection[j].A;
				++vx_count;
				break;
			case 1: 		 				// y-axis connection
				vy += connection[j].Flux(phase,cell,liquid,g)/connection[j].A;
				++vy_count;				
				break;
			default: 							// z-axis connection
				assert(connection[j].dir == 2);
				vz += connection[j].Flux(phase,cell,liquid,g)/connection[j].A;
				++vz_count;				
				break;
			}
		}
	vx /= max(vx_count,1);
	vy /= max(vy_count,1);
	vz /= max(vz_count,1);
	return tuple(vx,vy,vz);	
	}

void WriteVelocities(string file_name,Cell[] cell,Connection[] connection,Liquid[] liquid,double g) {
	// write results to output file
	double vx_w,vy_w,vz_w,vx_n,vy_n,vz_n;
	auto output_file = File(file_name, "w");
	output_file.writeln("cell","\t","x","\t","y","\t","z","\t","vx_w","\t","vy_w","\t","vz_w","\t","vx_napl","\t","vy_napl","\t","vz_napl");
	for (int i = 0; i < cell.length; ++i) {
		auto v_w = Velocity(i,0,connection,cell,liquid,g);
		vx_w = v_w[0];
		vy_w = v_w[1];
		vz_w = v_w[2];
		auto v_n = Velocity(i,1,connection,cell,liquid,g);
		vx_n = v_n[0];
		vy_n = v_n[1];
		vz_n = v_n[2];
		output_file.writeln(i,"\t",cell[i].x,"\t",cell[i].y,"\t",cell[i].z,"\t",vx_w,"\t",vy_w,"\t",vz_w,"\t",vx_n,"\t",vy_n,"\t",vz_n);
		}
	output_file.close();
	}	
	
void WriteFluxes(Connection[] connection,Cell[] cell,Liquid[] liquid,double g) {
	// write initial fluxes to file (to check model setup)
	double flux_w,flux_n;
	auto output_file = File("init_fluxes.txt", "w");
	output_file.writeln("connection","\t","cell_1","\t","cell_2","\t","Qw","\t","Qn");
	for (int i = 0; i < connection.length; ++i) {
		flux_w = connection[i].Flux(0,cell,liquid,g);
		flux_n = connection[i].Flux(1,cell,liquid,g);		
		output_file.writeln(i,"\t",connection[i].icell_0,"\t",connection[i].icell_1,"\t",flux_w,"\t",flux_n);
		}
	output_file.close();
	}	
	
void main(){

	double time = 0.0;
	Grid grid;
	Model model;
	Source source_w,source_n;
	Liquid[] liquid;
	Cell[] cell;
	Connection[] connection;
	const double g = 9.807;				// gravitational acceleration
	
	// read in liquid properties
	liquid = LiquidProps();
	writeln("Read liquid properties.");	
	
	// read in grid geometry
	grid = new Grid();
	writeln("Read grid geometry.");	
	
	// read in/set up cell set
	cell = ConstructCells(grid,liquid,g);
	for (int i = 0; i < cell.length; ++i){cell[i].P = cell[i].Pressure(liquid,g);}
	writeln("Set up cell set.");
	
	// establish cell connections
	connection = Connect(grid);
	cell = CellIndices(cell,connection);
	writeln("Set up cell connections.");
	
	// read in source terms
	source_w = new Source(cell.length,"source_water.txt");
	source_n = new Source(cell.length,"source_NAPL.txt");	
	writeln("Set up source terms.");	

	// set up model object
	model = new Model(cell.length,g);
	writeln("Set up simulation constraints.");	
	
	// run full simulation and generate output files
	writeln("Running simulation ...");		
	model.Simulate(cell,connection,liquid,source_w,source_n,time);
	
	writeln("Finished.");	
	
}