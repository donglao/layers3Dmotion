#include "ConnectedComponentDetection.h"


int ConnectedComponentDetection::allocate( int XSize, int YSize )
{
	this->XSize=XSize;
	this->YSize=YSize;
	GridSize = XSize*YSize;

	if ( !( ccLabel = new int[GridSize] ) || !( stack = new int[GridSize] ) ) {
		deallocate();
		return 0;
	}

	return 1;
}


void ConnectedComponentDetection::deallocate( )
{
	delete[] ccLabel; ccLabel=0;
	delete[] stack; stack=0;

	return;
}


/*
 * finds and labels cc of the pixels s.t. label==1
 */
void ConnectedComponentDetection::findConnectedComponents( const char *label )
{
	for (int p=0; p<GridSize; p++)  {
		if (label[p]!=0) ccLabel[p]=-1;
		else ccLabel[p]=0;
	}

	Ncomponents=0;

	for (int p=0; p<GridSize; p++) {
		if (ccLabel[p] == -1) {
			Ncomponents++;
			fill( p, Ncomponents, label );
		}
	}

	return;
}


void ConnectedComponentDetection::fill( int seed, int id, const char *label )
{
	int p, x, y;
	int seed_label = label[seed];
	int *top_stack=stack;

#define push(p) \
	if ( ccLabel[p] == -1 && label[p] == seed_label ) { \
		*top_stack++=(p);                               \
		ccLabel[p]=id;                                  \
	}

#define pop() ( top_stack!=stack ?  (*--top_stack) : -1 )

	push( seed );
	while ( (p=pop())!=-1 ) {
		x=p%XSize; y=p/XSize;
		if ( x < XSize-1 ) { push(p+1); }
		if ( x > 0 )       { push(p-1); }
		if ( y < YSize-1 ) { push(p+XSize); }
		if ( y > 0 )       { push(p-XSize); }
	}
	
#undef push
#undef pop

	return;
}