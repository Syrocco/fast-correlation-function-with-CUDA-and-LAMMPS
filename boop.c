/*
   Bond orientational order parameter calculation.
 */

#include<math.h>
#include<stdbool.h>
#include<getopt.h>
#include<string.h>
#include"parser.h"

#define PATH_SIZE 256
#define NMAX 10      //Maximum number of boops that can be computed in one call
#define NTYPEMAX 10  //Maximum number of allowed particle types 

bool vflag=0;  //Verbose flag (set by -v option)
bool qflag=0;  //Quiet flag (set by -q option)

struct List{
	int val;
	struct List* next;
};

int type_index(int type,int* allowed_types,size_t ntypes){
	/*
	 * Returns index of type in allowed_types.
	 * Return -1 if type is not found in the ntypes first values.
	 */
	for(unsigned int k=0;k<ntypes;++k){
		if(type==allowed_types[k]){
			return k;
		}
	}
	return -1;
}

double max_cutoff(double (*cutoff)[NTYPEMAX],size_t ntypes){
	double m=0;
	for(size_t i=0;i<ntypes;++i){
		for(size_t j=i;j<ntypes;++j){
			if(cutoff[i][j]>m){m=cutoff[i][j];}
		}
	}
	return sqrt(m);
}

void boop_framecl(int* n,int nn,int* allowed_types,int ntypes,double cutoff[][ntypes],bool savecomplex,Dump* dump,Dump* dumpout){
	/*
	 * Computes bond order parameters for qi for specified i's, only considering particles of specified types. Adds the qi property to atoms and stores the updated frame in fileout.
	 * n contains the nn values of i for qi computations.
	 * allowed_types specifies the ntypes types for which boop is to be computed.
	 * cutoff is a ntype*ntypes matrix containing the square of the distance below which two particles of types i and j or considered neighbours.
	 * Uses cell list to speed up calculation : should not be used on small systems.
	 */
	if(vflag){printf("Frame: %d\n",dump->curframe);}
	int natoms=get_natoms(dump);


	double (*qR)[natoms] = calloc(nn, sizeof(*qR));
	double (*qI)[natoms] = calloc(nn, sizeof(*qI));
	double (*qn)[natoms] = calloc(nn, sizeof(*qn));
	double (*qA)[natoms] = calloc(nn, sizeof(*qA));
	int *nb = calloc(natoms, sizeof(int));
	double* x = calloc(natoms, sizeof(double));
	double* y = calloc(natoms, sizeof(double));
	int* type = calloc(natoms, sizeof(int));
	
	get_doubleatomprop("x",x,natoms,dump);
	get_doubleatomprop("y",y,natoms,dump);
	get_intatomprop("type",type,natoms,dump);
	char hboxx;
	double lboxx=get_boxx(&hboxx,1,dump);
	char hboxy;
	double lboxy=get_boxy(&hboxy,1,dump);
	//Build cell list
	double cmax=max_cutoff(cutoff,ntypes);
	int ncx=lboxx/cmax;   //Number of cells
	int ncy=lboxy/cmax;
	if(ncx<3||ncy<3){
		printf("[boop_framecl]: System is too small for cell list optimisation. Use -L option to disable it.\n");
		exit(1);
	}
	double lcx=lboxx/ncx; //Length of cells
	double lcy=lboxy/ncy;
	struct List* cell_list[ncx][ncy];
	for(int i=0;i<ncx;++i){
		for(int j=0;j<ncy;++j){
			cell_list[i][j]=NULL;
		}
	}
	//Fill cell list
	struct List* cell = calloc(natoms, sizeof(struct List));
	for(int i=0;i<natoms;++i){
		if(type_index(type[i],allowed_types,ntypes)>=0){
			int indx=x[i]/lcx;
			int indy=y[i]/lcy;
			cell[i].val=i;
			cell[i].next=cell_list[indx][indy];
			cell_list[indx][indy]=&cell[i];
		}
	}
	//Compute boops using cell list
	for(int i=0;i<natoms;++i){
		int indti=type_index(type[i],allowed_types,ntypes);
		if(indti<0){continue;}
		int cellx=x[i]/lcx;
		int celly=y[i]/lcy;
		for(int a=-1;a<=1;++a){
			int indx=(cellx+a+ncx)%ncx;
			//Don't wrapp cells if boundaries are not periodic
			if(hboxx!='p' && cellx+a!=indx){continue;}
			for(int b=-1;b<=1;++b){
				int indy=(celly+b+ncy)%ncy;
				if(hboxy!='p' && celly+b!=indy){continue;}
				struct List* ll=cell_list[indx][indy];
				while(ll){
					int j=ll->val;
					int indtj=type_index(type[j],allowed_types,ntypes);
					if(j<=i){ll=ll->next;continue;}
					double dx=x[j]-x[i];
					if(hboxx=='p'){
						if(dx<-lboxx/2){dx+=lboxx;}
						if(dx>lboxx/2){dx-=lboxx;}
					}
					double dy=y[j]-y[i];
					if(hboxy=='p'){
						if(dy<-lboxy/2){dy+=lboxy;}
						if(dy>lboxy/2){dy-=lboxy;}
					}
					double r2=dx*dx+dy*dy;
					if(r2<cutoff[indti][indtj]){
						++nb[i];
						++nb[j];
						double theta=atan2(dy,dx);
						for(int k=0;k<nn;++k){
							double cc=cos(n[k]*theta);
							double ss=sin(n[k]*theta);
							int fact=n[k]%2?-1:1;
							qR[k][i]+=cc;
							qI[k][i]+=ss;
							qR[k][j]+=cc*fact;
							qI[k][j]+=ss*fact;
						}
					}
					ll=ll->next;
				}
			}
		}
		for(int k=0;k<nn;++k){
			if(nb[i]==0){
				qn[k][i]=0;
				continue;
			}
			qR[k][i]/=nb[i];
			qI[k][i]/=nb[i];
			qn[k][i]=sqrt(qR[k][i]*qR[k][i]+qI[k][i]*qI[k][i]);
			qA[k][i] = atan2(qI[k][i], qR[k][i]);
		}
	}
	add_intatomprop(nb,"nb",dumpout,dump);
	for(int k=0;k<nn;++k){
		char prop[]="q0";
		prop[1]+=n[k];
		add_doubleatomprop(qn[k],prop,dumpout,dump);
		if (savecomplex)
		{
			char prop2[] = "qR0";
			prop2[2] += n[k];
			add_doubleatomprop(qR[k],prop2,dumpout,dump);
			prop2[1] = 'I';
			add_doubleatomprop(qI[k],prop2,dumpout,dump);
			prop2[1] = 'A';
			add_doubleatomprop(qA[k],prop2,dumpout,dump);
		}
	}
	write_frame(dumpout);
	free(qR);
	free(qI);
	free(qn);
	free(qA);
	free(nb);
	free(x);
	free(y);
	free(type);
	free(cell);
}

void boop_frame(int* n,int nn,int* allowed_types,int ntypes,double cutoff[][ntypes],bool savecomplex,Dump* dump,Dump* dumpout){
	/*
	 * Computes bond order parameters for qi for specified i's, only considering particles of specified types. Adds the qi property to atoms and stores the updated frame in fileout.
	 * n contains the nn values of i for qi computations.
	 * allowed_types specifies the ntypes types for which boop is to be computed.
	 * cutoff is a ntype*ntypes matrix containing the square of the distance below which two particles of types i and j or considered neighbours.
	 */
	if(vflag){printf("Frame: %d\n",dump->curframe);}
	int natoms=get_natoms(dump);

	double (*qR)[natoms] = calloc(nn, sizeof(*qR));
	double (*qI)[natoms] = calloc(nn, sizeof(*qI));
	double (*qn)[natoms] = calloc(nn, sizeof(*qn));
	
	int *nb = calloc(natoms, sizeof(int));
	double* x = calloc(natoms, sizeof(double));
	double* y = calloc(natoms, sizeof(double));
	int* type = calloc(natoms, sizeof(int));

	get_doubleatomprop("x",x,natoms,dump);
	get_doubleatomprop("y",y,natoms,dump);
	get_intatomprop("type",type,natoms,dump);
	char hboxx;
	double lboxx=get_boxx(&hboxx,1,dump);
	char hboxy;
	double lboxy=get_boxy(&hboxy,1,dump);
	for(int i=0;i<natoms;++i){
		int indti=type_index(type[i],allowed_types,ntypes);
		if(indti<0){continue;}
		for(int j=i+1;j<natoms;++j){
			int indtj=type_index(type[j],allowed_types,ntypes);
			if(indtj<0){continue;}
			double dx=x[j]-x[i];
			if(hboxx=='p'){
				if(dx<-lboxx/2){dx+=lboxx;}
				if(dx>lboxx/2){dx-=lboxx;}
			}
			double dy=y[j]-y[i];
			if(hboxy=='p'){
				if(dy<-lboxy/2){dy+=lboxy;}
				if(dy>lboxy/2){dy-=lboxy;}
			}
			double r2=dx*dx+dy*dy;
			if(r2<cutoff[indti][indtj]){
				++nb[i];
				++nb[j];
				double theta=atan2(dy,dx);
				for(int k=0;k<nn;++k){
					double cc=cos(n[k]*theta);
					double ss=sin(n[k]*theta);
					int fact=n[k]%2?-1:1;
					qR[k][i]+=cc;
					qI[k][i]+=ss;
					qR[k][j]+=cc*fact;
					qI[k][j]+=ss*fact;
				}
			}
		}
		for(int k=0;k<nn;++k){
			qR[k][i]/=nb[i];
			qI[k][i]/=nb[i];
			qn[k][i]=sqrt(qR[k][i]*qR[k][i]+qI[k][i]*qI[k][i]);
		}
	}
	add_intatomprop(nb,"nb",dumpout,dump);
	for(int k=0;k<nn;++k){
		char prop[]="q0";
		prop[1]+=n[k];
		add_doubleatomprop(qn[k],prop,dumpout,dump);
		if (savecomplex)
		{
			char prop2[] = "qR0";
			prop2[2] += n[k];
			add_doubleatomprop(qR[k],prop2,dumpout,dump);
			prop2[1] = 'I';
			add_doubleatomprop(qI[k],prop2,dumpout,dump);
		}
	 }
	write_frame(dumpout);
	free(qR);
	free(qI);
	free(qn);
	free(nb);
	free(x);
	free(y);
	free(type);
}


int main(int argc,char** argv){
	int n[NMAX];
	int nn=0;
	int ntypes=0;
	int allowed_types[NTYPEMAX];
	double cutoff[NTYPEMAX][NTYPEMAX];
	bool c_encountered=0; //Set after -c: no more -t accepted
	bool Lflag=0;    //If set by -L, no cell list optimisation
	bool savecomplex = 0; // If set by --save-complex, save real and imaginary part of boops in addition to amplitude.
	//Parse commandline arguments
	struct option longopt[]={
		{"type",required_argument,NULL,'t'},
		{"cutoffs",required_argument,NULL,'c'},
		{"cell-list-off",no_argument,NULL,'L'},
		{"save-complex", no_argument, NULL, 1000},
		{"verbose",no_argument,NULL,'v'},
		{"quiet",no_argument,NULL,'q'},
		{"help",no_argument,NULL,'h'},
		{0,0,0,0}};
	int c;
	while((c=getopt_long(argc,argv,"n:t:c:Lvqh",longopt,NULL))!=-1){
		switch(c){
			case 'n':
				if(nn<NMAX){sscanf(optarg,"%d",&n[nn]);}
				else{
					printf("[boop]: A maximum of %d boops can be computed in one call. Ignoring extra values.\n",NMAX);
					break;
				}
				++nn;
				break;
			case 't':
				if(c_encountered){
					printf("-t --type options encountered after -c --cutoffs are ignored\n");
					break;
				}
				if(ntypes<NTYPEMAX){sscanf(optarg,"%d",&allowed_types[ntypes]);}
				else{
					printf("[boop]: A maximum of %d particle types can be handled in a single call. Ignoring extra values.\n",NTYPEMAX);
					break;
				}
				++ntypes;
				break;
			case 'c':
				; //Empty statement following label
				char* token=strtok(optarg," ");
				while(token!=NULL){
					int typei;
					int typej;
					double cut;
					sscanf(token,"%d:%d:%lf",&typei,&typej,&cut);
					int indti=type_index(typei,allowed_types,ntypes);
					int indtj=type_index(typej,allowed_types,ntypes);
					if(indti>=0 && indtj>=0){
						cutoff[indti][indtj]=cut*cut;
						cutoff[indtj][indti]=cut*cut;
					}
					else{
						printf("Unknown type in '%s'. Exit !\n",token);
						exit(1);
					}
					token=strtok(NULL," ");
				}
				c_encountered=1;
				break;
			case 'L':
				Lflag=1;
				break;
			case 1000:
				savecomplex = 1;
				break;
			case 'v':
				vflag=1;
				break;
			case 'q':
				qflag=1;
				break;
			case 'h':
				printf(
 " * Usage: ./boop [OPTIONS] SOURCE DESTINATION\n * Description: computes boop of specified order using position of atoms in LAMMPS-formatted dump SOURCE and saves the results in a LAMMPS-formatted dump DESTINATION. SOURCE and DESTINATION must be different.\n * Options:\n * -n [int]: add a value to the list of boops to be computed\n * -t --type [int]: add a value to the list of particles types to be considered\n * -c --cutoffs[string]: specifies cutoff distances. Format: \"typei:typej:cutoff typek:typel:cutoff ...\". Types must have been previously declared with -t --type option. Types provided after the -c --cutoff option are ignored.\n * -L --cell-list-off: disable cell list optimisation (necessary for small systems).\n * --save-complex: save real and imaginary part of boops in addition to amplitude.\n * -v: verbose.\n * -q: quiet.\n * -h: display this message and exit.\n"
				 );
				exit(0);
		}
	}
	if(optind>argc-2){
		printf("Usage: ./boop [OPTIONS] SOURCE DESTINATION\n");
		exit(0);
	}
	if(nn==0){
		printf("No boop calculation requested. Use -n option to specify desired order (may be used multiple times).\n");
		exit(0);
	}
	if(ntypes==0){
		printf("No allowed types specified. Use -t option to specify the types of particles to be considered (may be used multiple times).\n");
		exit(0);
	}
	char* path_in=argv[optind];
	char* path_out=argv[optind+1];
	Dump* dump=dump_open(path_in,'r');
	Dump* dumpout=dump_open(path_out,'w');
	unsigned int nframes=dump->nframes;
	if(!qflag){
		printf("Bond order parameters to be computed: ");
		for(int i=0;i<nn;++i){printf("q%d ",n[i]);}
		printf("\n");
		printf("Allowed particle types: ");
		for(int i=0;i<ntypes;++i){printf("%d ",allowed_types[i]);}
		printf("\n");
		printf("Cutoffs:\n");
		for(int i=0;i<ntypes;++i){
			for(int j=i;j<ntypes;++j){
				printf("%d-%d: %lf\n",allowed_types[i],allowed_types[j],sqrt(cutoff[i][j]));
			}
		}
		if(Lflag)
			printf("Cell-list optimisation disabled.\n");
		if (savecomplex)
			printf("Saving real and imaginary part of boops.\n");
		printf("Input: %s\nNumber of frames: %d\n",path_in,nframes);
		printf("Output: %s\n",path_out);
		printf("\nComputing boops...\n");
	}
	if(Lflag){
		for(unsigned int f=0;f<nframes;++f){
			jump_to_frame(f,dump);
			boop_frame(n,nn,allowed_types,ntypes,cutoff,savecomplex,dump,dumpout);
			
		}
	}
	else{
		//Compute boops for each frame
		for(unsigned int f=0;f<nframes;++f){
			jump_to_frame(f,dump);
			boop_framecl(n,nn,allowed_types,ntypes,cutoff,savecomplex,dump,dumpout);
			
		}
	}
	dump_close(dumpout);
	dump_close(dump);
}
