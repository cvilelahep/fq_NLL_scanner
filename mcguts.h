/*****************************************************************

$Id: mcguts.h,v 1.1 2012/06/16 15:39:11 shimpeit Exp $

$Log: mcguts.h,v $
Revision 1.1  2012/06/16 15:39:11  shimpeit
Importing atmpd/branches/dev/gfortran w/ fiTQun

Revision 1.1  2012/03/27 00:50:06  wilking
Add first version of pi0 fitter.

Revision 1.1  2000/03/07 20:15:02  habig
moved mcguts to libmcutil.a


mcguts.h: Header file for mcguts.c
-Matt Earl : August 15, 1997

******************************************************************/

enum what_it_is {INCOMING = 0, OUTGOING, SECONDARY};
enum fs {PRIMARYINT = -1, UNKNOWN, DECAY, ESCAPE, ABSORB, CHARGEX, STOP, EMSHOWER, HADINT};

#ifndef NMAXTRACKS
#define NMAXTRACKS 1000
#endif


struct Mc_particle
{
  int ntracks;                       /* Number of particles */
  int nincoming;                     /* Number of incoming particles */
  int noutgoing;                     /* Number of outgoing particles */
  int nsecondary;                    /* Number of secondary particles */
  int event_type;
  float primary_vertex[3];           /* Primary Vertex */
  float primary_dir[3];              /* Primary Direction */
  int type[NMAXTRACKS];              /* PDG code for the particle */
  int charge[NMAXTRACKS]; 
  int origin[NMAXTRACKS];            /* incoming, outgoing, or secondary */
  int final_state[NMAXTRACKS];       /* from fs above */
  int tracked[NMAXTRACKS];           /* particle was tracked in detector sim */
  int parent[NMAXTRACKS];
  float mass[NMAXTRACKS];            /* mass of particle in MeV/c^2 */
  float momentum[NMAXTRACKS];        /* momentum of particle in MeV/c */
  float energy[NMAXTRACKS];          /* energy of particle in MeV */
  float beta[NMAXTRACKS];
  float dir[NMAXTRACKS][3];
  float initial_vertex[NMAXTRACKS][3];
  float final_vertex[NMAXTRACKS][3];
  float time[NMAXTRACKS];            /* time particle created */
  int   ancestor[NMAXTRACKS];        /* the primary vector which eventually
					spawned this particle */
};
extern struct Mc_particle mc_particle;

/*void mcguts(int *ntracks);*/
char *ParticleName(int pdg_code);
int charge(int pdg_code);
char *EventType(int event_code);
const char *InteractionType(void);
void PrintMCInfo(void);
void PrintMCTrackInfo(int);
void SetmcgutsNEUT(void);
void SetmcgutsNUANCE(void);
void SetmcgutsNoVerbose(void);
void SetmcgutsVerbose(void);
