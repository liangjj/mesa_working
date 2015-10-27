/* Take a molecule structure and print it to a file.  Do it by violating
   all the rules and sticking our nose the stinking low-level data
   structure! */ 

#include <stdio.h>
#include <avs/avs.h>
#define CHEM_APPL
#include <CHEMmol.h>
#include <CHEMatom.h>
#include <CHEMquant.h>
#include <mol_udat.h>
#include <atom_udat.h>

#define AtomNameMask (1<<0)
#define AtomInumberMask (1<<1)
#define AtomColorMask (1<<2)
#define AtomRadiusMask (1<<3)
#define AtomCandBMask (1<<4)
#define AtomXyzMask (1<<5)
#define AtomHybridizationMask (1<<1)
#define AtomAtomicNumberMask (1<<2)
#define AtomChargeMask (1<<6)
#define TRUE 1
#define FALSE 0
#define ERROR -1

void AVSinit_modules()
{
  int voyeur();

  AVSmodule_from_desc(voyeur);

}

int voyeur()
{
  int peep_through_window(),param=0;
  char *directory=NULL,getcwd(),dir[256];

  AVSset_module_name("Peeping Tom",MODULE_RENDER);
  AVScreate_input_port("molecule","molecule",REQUIRED);

  directory=getcwd(NULL,256);
  sprintf(dir,"%s\n",directory);

  param=AVSadd_parameter("Output File","string",dir,NULL,".DMP");
  AVSconnect_widget(param,"browser");

  AVSset_compute_proc(peep_through_window);
}

int peep_through_window(inmol,filename)
CHEMmolecule *inmol;
char *filename;
{
   /*  What we're gonna do is just like CHEM_molecule_print followed by a 
       formatted dump of the molecule's user_data.  */

  FILE *fp;
  char *pos_dot=NULL;
  int dummy=0;
  int i;
  int struct_dims[NUMCHEMTYPES];
  CHEMatom *theatom;
  CHEMquantum *thequant;

  if (filename == NULL) return(FALSE);
  if ((pos_dot=strrchr(filename,'.'))==NULL) return(FALSE);
  if (strcmp(pos_dot,".DMP")) 
    {
      AVSwarning("Filename has no .DMP extension");
      return(FALSE);
    }
  if (CHEMgen_util_update_molecule(inmol,SINGLE,struct_dims,&dummy,&dummy))
    AVSerror("update molecule");

  if ((fp=fopen(filename,"w"))==NULL)
    {
      AVSwarning("Unable to open output file");
      return (FALSE);
    }

  /* here's where we cheat, because I don't want to use the accessor
     functions... blech */

  fprintf(fp,"CHEMmolecule:");
  fprintf(fp,"%s:",inmol->name);
  fprintf(fp,"%s:",((inmol->units)==0)?"angstroms":"bohr");
  fprintf(fp,"%d:%d:%d.\n",(inmol->natom),(inmol->nc_unit),(inmol->nquant));
  dump_mol_user_data(fp,inmol);
  

  theatom=inmol->atom;
  for (i=0;i<(inmol->natom)&&theatom != NULL;i++)
    {
      sniff_atom(fp,theatom);
      dump_atom_user_data(fp,theatom);
      theatom=theatom->next;
    }
  

  thequant=inmol->quant;
  for (i=0;i<inmol->nquant && thequant!= NULL;i++)
    {
      sniff_quant(fp,thequant);
      thequant=thequant->next;
    }
  fclose(fp);
  return (TRUE);
}

sniff_atom(fp,atom)
CHEMatom *atom;
FILE *fp;
{
  fprintf(fp,"  CHEMatom:");
  fprintf(fp,"%d:%s:%u:%lf,%lf,%lf:%f\n",
	  atom->index_number,atom->name,atom->color,atom->x,atom->y,atom->z,
	  atom->radius);
}

dump_atom_user_data(fp,atom)
FILE *fp;
CHEMatom *atom;
{

  MSIatom_udat *ud;

  if (atom->user_data == NULL) 
    {
      fprintf(fp,"    MSIatom_udat not present\n");
      return;
    }
  ud = atom->user_data;
  fprintf(fp,"     MSIatom_udat:");
  fprintf(fp,"%s:%s:",ud->version,ud->symbol);
  fprintf(fp,"%d:%d:%f:%d:%d:%d:%f:%d:%d\n",ud->setatominfo,
	  ud->atomic_number,ud->atomic_weight,ud->hybridization,ud->reptype,
	  ud->picked_color,ud->charge,ud->number_of_nubs,ud->nub_color);

}


dump_mol_user_data(fp,mol)
CHEMmolecule *mol;
FILE *fp;
{
  MSImolecule_udat *ud;

  if (mol->user_data == NULL)
    {
      fprintf(fp,"   MSImolecule_udat not present\n");
      return ; 
    }
  ud=mol->user_data;

  fprintf(fp,"   MSImolecule_udat:");
  fprintf(fp,"%s:%s:",ud->version,ud->creator);
  fprintf(fp,"%d:%d:%d:%d:%f,%f,%f\n",ud->chem_molatominfo,ud->msi_molatominfo,
	  ud->chem_molinfo,ud->msi_molinfo,ud->dipole[0],ud->dipole[1],
	  ud->dipole[2]);
}

sniff_quant(fp,q)
CHEMquantum *q;
FILE *fp;
{

  CHEMq_page *qp;
  int i;

  fprintf(fp,"  CHEMquantum:");

  fprintf(fp,"%s:",q->name);
  fprintf(fp,"%d:%d:%d:%d:%d:%d:%d:%d:%d:%d:%d\n",
	  q->nbasis,q->ich,q->mul,q->ne,q->na,q->nb,q->scftype,
	  q->corrtype,q->npage,q->nshell,q->ngauss);

  qp=q->pl;
  for (i=0;i<q->npage && qp!=NULL;i++)
    {
      fprintf(fp,"    CHEMq_page:");
      fprintf(fp,"%s:%d:%d\n",qp->name,qp->mbasis,qp->nbasis);
      qp=qp->next;
    }
}
