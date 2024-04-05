/*
 * Usage: ./listneigh [OPTIONS] DUMP OUTPUT
 * Description: List neighbours of particles of specified type in specified frames of LAMMPS-formatted DUMP and write neighbours lists into OUTPUT. OUTPUT contains the selected frames, with a new ITEM: NEIGHBOURS.
 * Options:
 * -f --frame[int]: Specify one frame for which neighbours are to be determined. Can be used multiple times (at most MAXFRAME). To select a range of frame, prefer -r --range option.
 * -r --range[string]: Specify a range of frames. The string must be "start:step:stop" or "*" in which case all frames are selected. If stop is set to -1, it is replaced by the number of frames in the dump file.
 * -t --type[int]: Specify one particle type to be considered in the computation. Can be used multiple time (at most NTYPEMAX).
 * -c --cutoffs[string]: Specify the cutoffs for each type pair of selected types. The format of the string must be "type0:type1:cutoff01 ... typei:typej:cutoffij"
 * -L --cell-list_off: disable cell list optimisation (necessary for small systems).
 * -v --verbose: verbose.
 * -q --quiet: quiet.
 * -h --help: display this message and exit.
 */

#include<stdio.h>
#include<math.h>
#include<getopt.h>
#include"parser.h"

#define NFRAMEMAX 10  //Maximum number of frames that can be specified with -f
#define NTYPEMAX 10   //Maximum number of allowed types that can be specified wth -t
#define NEIGHMAX 15   //Maximum number of neighbours for one particle

bool vflag=0;	//Verbose flag, set by -v
bool qflag=0;	//Quiet flag, set by -q

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

void write_neigh(int neigh[][NEIGHMAX],double cutoffs[][NTYPEMAX],int natoms,int nb[natoms],int type[natoms],int allowed_types[NTYPEMAX],int ntypes,Dump* dump,Dump* dumpout){
	/*
	 * Copies the current frame of dump in dumpout.
	 * Adds nb property to atoms (number of neighbours).
	 * Adds an ITEM: NEIGHBOURS section with neighbours lists.
	 * ITEM: NEIGHBOURS type0:type1:cutoff ... typei:typej:cutoff
	 * 0: neigh0 neigh1 ... neighk
	 * ...
	 */
	add_intatomprop(nb,"nb",dumpout,dump);
	// Forge header
	int headsize=255;
	char* headbuff=(char*)malloc((headsize+1)*sizeof(char));
	if(headbuff==NULL){
		printf("[write_neigh]: Failed to allocate neigh buffer of size %d.\n",headsize+1);
		exit(1);
	}
	int charcount=0;
	int writtenchar=0;
	writtenchar=snprintf(headbuff+charcount,headsize+1-charcount,"ITEM: NEIGHBOURS ");
	// Grow buffer if previous write filled it, and write again
	while(writtenchar+charcount>=headsize){
		headsize*=2;
		headbuff=(char*)realloc(headbuff,(headsize+1)*sizeof(char));
		if(headbuff==NULL){
			printf("[write_neigh]: Failed to allocate neigh buffer of size %d.\n",headsize+1);
			exit(1);
		}
		writtenchar=snprintf(headbuff+charcount,headsize+1-charcount,"ITEM: NEIGHBOURS ");
	}
	charcount+=writtenchar;
	for(int i=0;i<ntypes;++i){
		for(int j=i;j<ntypes;++j){
			writtenchar=snprintf(headbuff+charcount,headsize+1-charcount,"%d:%d:%lf ",allowed_types[i],allowed_types[j],sqrt(cutoffs[i][j]));
			while(writtenchar+charcount>=headsize){
				headsize*=2;
				headbuff=(char*)realloc(headbuff,(headsize+1)*sizeof(char));
				if(headbuff==NULL){
					printf("[write_neigh]: Failed to allocate neigh buffer of size %d.\n",headsize+1);
					exit(1);
				}
				writtenchar=snprintf(headbuff+charcount,headsize+1-charcount,"%d:%d:%lf ",allowed_types[i],allowed_types[j],sqrt(cutoffs[i][j]));
			}
			charcount+=writtenchar;
		}
	}
	headbuff[charcount-1]='\n'; // Replace last space with newline
	headsize=charcount;
	// Forge section content
	charcount=0;
	int neighsize=dumpout->framesize;
	char* buffneigh=(char*)malloc((neighsize+1)*sizeof(char));
	if(buffneigh==NULL){
		printf("[write_neigh]: Failed to allocate neigh buffer of size %d.\n",neighsize+1);
		exit(1);
	}
	for(int i=0;i<natoms;++i){
		if(type_index(type[i],allowed_types,ntypes)!=-1){
			// Write id of the considered atom
			writtenchar=snprintf(buffneigh+charcount,neighsize+1-charcount,"%d: ",i);
			while(charcount+writtenchar>=neighsize){
				neighsize*=2;
				buffneigh=(char*)realloc(buffneigh,(neighsize+1)*sizeof(char));
				if(buffneigh==NULL){
					printf("[write_neigh]: Failed to allocate neigh buffer of sier %d.\n",neighsize);
					exit(1);
				}
				writtenchar=snprintf(buffneigh+charcount,neighsize+1-charcount,"%d: ",i);
			}
			charcount+=writtenchar;
			// Add its neighbours to the line
			for(int k=0;k<nb[i];++k){
				writtenchar=snprintf(buffneigh+charcount,neighsize+1-charcount,"%d ",neigh[i][k]);
				while(charcount+writtenchar>=neighsize){
					neighsize*=2;
					buffneigh=(char*)realloc(buffneigh,(neighsize+1)*sizeof(char));
					if(buffneigh==NULL){
						printf("[write_neigh]: Failed to allocate neigh buffer of sier %d.\n",neighsize);
						exit(1);
					}
					writtenchar=snprintf(buffneigh+charcount,neighsize+1-charcount,"%d ",neigh[i][k]);
				}
				charcount+=writtenchar;
			}
			buffneigh[charcount-1]='\n';
		}
	}
	neighsize=charcount;
	// Write the new section to dumpout
	add_section(headbuff,headsize,buffneigh,neighsize,dumpout,dump);
	free(headbuff);
	free(buffneigh);
	write_frame(dumpout);
}

void listneigh_frame(int* allowed_types,int ntypes,double cutoff[][ntypes],Dump* dump,Dump* dumpout){
	/*
	 * Run when -L option is used.
	 * Computes neighbours list for the current frame of dump by listing all pairs of particles. For large systems, prefer listneigh_frame_cl that uses cell list optimisation.
	 * Writes the output to dumpout.
	 */
	if(vflag==1){printf("  Frame %d\n",dump->curframe);}
	int natoms=get_natoms(dump);
	int nb[natoms];   //Number of neighbours
	memset(nb,0,sizeof(nb));
	int neigh[natoms][NEIGHMAX];
	double x[natoms];
	double y[natoms];
	int type[natoms];
	get_doubleatomprop("x",x,natoms,dump);
	get_doubleatomprop("y",y,natoms,dump);
	get_intatomprop("type",type,natoms,dump);
	char hboxx_ptr[2];
	double lboxx=get_boxx(hboxx_ptr,2,dump);
	char hboxy_ptr[2];
	double lboxy=get_boxy(hboxy_ptr,2,dump);
	char hboxx=hboxx_ptr[0];
	char hboxy=hboxy_ptr[0];
	for(int i=0;i<natoms;++i){
		if(type_index(type[i],allowed_types,ntypes)>=0){
			for(int j=i+1;j<natoms;++j){
				if(type_index(type[j],allowed_types,ntypes)>=0){
					double dx=x[i]-x[j];
					double dy=y[i]-y[j];
					if(hboxx=='p'){
						if(dx<-lboxx/2){dx+=lboxx;}
						if(dx>lboxx/2){dx-=lboxx;}
					}
					if(hboxy=='p'){
						if(dy<-lboxy/2){dy+=lboxy;}
						if(dy>lboxy/2){dy-=lboxy;}
					}
					if(dx*dx+dy*dy<=cutoff[type[i]][type[j]]){
						if(nb[i]==NEIGHMAX){
							printf("[listneigh_frame]: particle %d has more than %d neihbours. Recompile this program with a larger NEIGHMAX.\n",i,NEIGHMAX);
							exit(1);
						}
						if(nb[j]==NEIGHMAX){
							printf("[listneigh_frame]: particle %d has more than %d neihbours. Recompile this program with a larger NEIGHMAX.\n",j,NEIGHMAX);
							exit(1);
						}
						neigh[i][nb[i]]=j;
						neigh[j][nb[j]]=i;
						++nb[i];
						++nb[j];
					}
				}
			}
		}
	}
	write_neigh(neigh,cutoff,natoms,nb,type,allowed_types,ntypes,dump,dumpout);
}

void listneigh_frame_cl(int* allowed_types,int ntypes,double cutoff[][ntypes],Dump* dump,Dump* dumpout){
	/*
	 * Computes neighbours list for the current frame of dump using cell lists. If the system size is too small to fit 3*3 cells, listneigh_frame must be used (use -L option).
	 * allowed_types contains the ntypes particles types to be considered in neighbour detection.
	 * i,j element of cutoff contains the square of the cutoff distance for particles of type allowed_type[i]-allowed_type[j].
	 * Writes the output to dumpout.
	 */
	if(vflag==1){printf("  Frame %d\n",dump->curframe);}
	int natoms=get_natoms(dump);
	int nb[natoms];   //Number of neighbours
	memset(nb,0,sizeof(nb));
	int neigh[natoms][NEIGHMAX];
	double x[natoms];
	double y[natoms];
	int type[natoms];
	get_doubleatomprop("x",x,natoms,dump);
	get_doubleatomprop("y",y,natoms,dump);
	get_intatomprop("type",type,natoms,dump);
	char hboxx_ptr[2];
	double lboxx=get_boxx(hboxx_ptr,2,dump);
	char hboxy_ptr[2];
	double lboxy=get_boxy(hboxy_ptr,2,dump);
	char hboxx=hboxx_ptr[0];
	char hboxy=hboxy_ptr[0];
	//Build cell list
	double cmax=max_cutoff(cutoff,ntypes);
	int ncx=lboxx/cmax;   //Number of cells
	int ncy=lboxy/cmax;
	if(ncx<3||ncy<3){
		printf("[listneigh_frame_cl]: System is too small for cell list optimisation. Use -L option to disable it.\n");
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
	struct List cell[natoms];
	for(int i=0;i<natoms;++i){
		if(type_index(type[i],allowed_types,ntypes)>=0){
			int indx=x[i]/lcx;
			int indy=y[i]/lcy;
			cell[i].val=i;
			cell[i].next=cell_list[indx][indy];
			cell_list[indx][indy]=&cell[i];
		}
	}
	//Find neighbours using cell list
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
						if(nb[i]==NEIGHMAX){
							printf("[listneigh_frame_cl]: particle %d has more than %d neihbours. Recompile this program with a larger NEIGHMAX.\n",i,NEIGHMAX);
							exit(1);
						}
						if(nb[j]==NEIGHMAX){
							printf("[listneigh_frame_cl]: particle %d has more than %d neihbours. Recompile this program with a larger NEIGHMAX.\n",j,NEIGHMAX);
							exit(1);
						}
						neigh[i][nb[i]]=j;
						neigh[j][nb[j]]=i;
						++nb[i];
						++nb[j];
					}
					ll=ll->next;
				}
			}
		}
	}
	write_neigh(neigh,cutoff,natoms,nb,type,allowed_types,ntypes,dump,dumpout);
}

int main(int argc, char** argv){
	int frames[NFRAMEMAX];
	int nframes=0;
	int allowed_types[NTYPEMAX];
	double cutoff[NTYPEMAX][NTYPEMAX];
	memset(cutoff,0,sizeof(cutoff));
	int ntypes=0;
	bool rflag=0;	//Set by -r, to use a range of frames
	bool selectall=0; //Set by -r *
	bool lflag=0;   //Set by -L, to disable cell list optimisation
	int rstart;
	int rstep;
	int rstop;

	struct option longopt[]={
		{"frame",required_argument,NULL,'f'},
		{"range",required_argument,NULL,'r'},
		{"type",required_argument,NULL,'t'},
		{"cutoffs",required_argument,NULL,'c'},
		{"cell-list-off",no_argument,NULL,'L'},
		{"verbose",no_argument,NULL,'v'},
		{"quiet",no_argument,NULL,'q'},
		{"help",no_argument,NULL,'h'},
		{0,0,0,0}};
	int c;
	while((c=getopt_long(argc,argv,"f:r:t:c:Lvqh",longopt,NULL))!=-1){
		switch(c){
			case 'f':
				if(rflag){
					printf("[listneigh]: -t --type is ignored when -r --range is used.\n");
					break;
				}
				sscanf(optarg,"%d",&frames[nframes]);
				++nframes;
				break;
			case 'r':
				if(rflag){printf("[listneigh]: Frame range already set. Overwriting.\n");}
				if(nframes>0){printf("[listneigh]: -t --type is ignored when -r --range is used.\n");}
				rflag=1;
				selectall=0;
				if(optarg[0]=='*'){
					selectall=1;
					break;
				}
				char* token=strtok(optarg,":");
				sscanf(token,"%d",&rstart);
				token=strtok(NULL,":");
				sscanf(token,"%d",&rstep);
				token=strtok(NULL,":");
				sscanf(token,"%d",&rstop);
				break;
			case 't':
				if(ntypes<NTYPEMAX){
					sscanf(optarg,"%d",&allowed_types[ntypes]);
					++ntypes;
				}
				else{
					printf("[listneigh]: A maximum of %d particle types can be considered in one call. Extra values are ignored.\n",NTYPEMAX);
				}
				break;
			case 'c':
				;
				int t0;
				int t1;
				double cut;
				char* tokenc=strtok(optarg,":");
				sscanf(tokenc,"%d",&t0);
				tokenc=strtok(NULL,":");
				sscanf(tokenc,"%d",&t1);
				tokenc=strtok(NULL,":");
				sscanf(tokenc,"%lf",&cut);
				if(t0>=NTYPEMAX){
					printf("[listneigh]: invalid type %d in -c option.\n",t0);
					exit(1);
				}
				if(t1>=NTYPEMAX){
					printf("[listneigh]: invalid type %d in -c option.\n",t1);
					exit(1);
				}
				cutoff[t0][t1]=cut*cut;
				cutoff[t1][t0]=cut*cut;
				break;
			case 'L':
				lflag=1;
				break;
			case 'v':
				vflag=1;
				break;
			case 'q':
				qflag=1;
				break;
			case 'h':
				printf(" * Usage: ./listneigh [OPTIONS] DUMP OUTPUT\n * Description: List neighbours of particles of specified type in specified frames of LAMMPS-formatted DUMP and write neighbours lists into OUTPUT. OUTPUT contains the selected frames, with a new ITEM: NEIGHBOURS.\n * Options:\n * -f --frame[int]: Specify one frame for which neighbours are to be determined. Can be used multiple times (at most MAXFRAME). To select a range of frame, prefer -r --range option.\n * -r --range[string]: Specify a range of frames. The string must be \"start:step:stop\" or \"*\", in which case all frames are selected. If stop is set to -1, it is replaced by the number of frames in the dump file.\n * -t --type[int]: Specify one particle type to be considered in the computation. Can be used multiple time (at most NTYPEMAX).\n * -c --cutoffs[string]: Specify the cutoffs for each type pair of selected types. The format of the string must be \"type0:type1:cutoff01 ... typei:typej:cutoffij\"\n * -L --cell-list-off: disable cell list optimisation (necessary for small systems).\n * -v --verbose: verbose.\n * -q --quiet: quiet.\n * -h --help: displays this message and exits.\n");
				exit(0);
		}
	}
	if(argc<optind+2){
		printf("[listneigh]: Usage: ./listneigh [OPTIONS] DUMP DUMPOUT\n");
		exit(0);
	}
	if(nframes==0 && rflag==0){
		printf("[listneigh]: No frame selected. Use -f --frame or -r --range option.\n");
		exit(0);
	}	
	if(ntypes==0){
		printf("[listneigh]: No particle type selected. Use -t --type option.\n");
		exit(0);
	}	
	Dump* dump=dump_open(argv[optind],'r');
	Dump* dumpout=dump_open(argv[optind+1],'w');
	if(selectall){
		rstart=0;
		rstep=1;
		rstop=dump->nframes;
	}
	if(rstop==-1){rstop=dump->nframes;}
	if(!qflag){
		printf("Input dump: %s\n  Found %d frames.\n",argv[optind],dump->nframes);
		printf("Output dump: %s\n",argv[optind+1]);
		printf("Selected frames: ");
		if(!rflag){for(int f=0;f<nframes;++f){printf("%d ",frames[f]);}}
		else{printf("range %d:%d:%d",rstart,rstep,rstop);}
		printf("\nParticle types considered: ");
		for(int i=0;i<ntypes;++i){printf("%d ",allowed_types[i]);}
		if(lflag==1){printf("\nCell list optimisation disabled.\n");}
		printf("\nComputing neighbours lists...\n");
	}
	if(!rflag){
		for(int f=0;f<nframes;++f){
			jump_to_frame(frames[f],dump);
			if(lflag==0){
				listneigh_frame_cl(allowed_types,ntypes,cutoff,dump,dumpout);
			}
			else{
				listneigh_frame(allowed_types,ntypes,cutoff,dump,dumpout);
			}
		}
	}
	else{
		for(int f=rstart;f<rstop;f+=rstep){
			jump_to_frame(f,dump);
			if(lflag==0){
				listneigh_frame_cl(allowed_types,ntypes,cutoff,dump,dumpout);
			}
			else{
				listneigh_frame(allowed_types,ntypes,cutoff,dump,dumpout);
			}
		}
	}
}

