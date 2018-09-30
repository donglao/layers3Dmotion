#ifndef CONNECTEDCOMPONENTDETECTION_H
#define CONNECTEDCOMPONENTDETECTION_H


class ConnectedComponentDetection {
	public:
		int Ncomponents;
		int *ccLabel;
		int XSize, YSize, GridSize;

	public:
		ConnectedComponentDetection() {
			Ncomponents=0;
			ccLabel=0;
			XSize=YSize=GridSize=0;
		}
		int  allocate(int XSize, int YSize);
		void deallocate();
		void findConnectedComponents( const char *label );
		void fill( int seed, int id, const char *label );

	private:
		int *stack;
};

#endif