/* This file is part of Cloudy and is copyright (C)1978-2025 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*SaveSpecial generate output for the save special command */
#include "cddefines.h"
#include "save.h"
#include "wind.h"
#include "opacity.h"
#include "dense.h"
#include "radius.h"
#include "colden.h"

/*SaveSpecial generate output for the save special command */
void SaveSpecial(FILE* ioPUN , 
  const char *chTime)
{
	/*long int i;*/

	DEBUG_ENTRY( "SaveSpecial()" );

	if( strncmp(chTime,"LAST",4) == 0 )
	{
		/* code to execute only after last zone */
		double wmean = 0.;
		if( colden.TotMassColl > 0. )
		{
			wmean = colden.wmas/SDIV(colden.TotMassColl);
		}

		fprintf(ioPUN,"# Final mean properties\n");
		fprintf(ioPUN,"MeanMolecularWeight\t");
		PrintE82(ioPUN , wmean);
		fprintf(ioPUN,"\n");

	}
	else
	{
		/* code to do for every zone */
		fprintf(ioPUN,"%.5e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n",
			radius.Radius ,
			wind.AccelCont ,
			wind.fmul ,
			opac.opacity_sct[1000],
			dense.eden , 
			dense.xMassDensity,
			dense.gas_phase[ipHYDROGEN] );
	}

	return;
}
