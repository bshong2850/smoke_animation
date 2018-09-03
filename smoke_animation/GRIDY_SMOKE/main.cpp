#include <iostream>
#include "../LIBRARY/GRIDY_SMOKE_3D.h"

#include "../LIBRARY/GRIDY_CONFIG.h"


int main(int argc, char **argv)
{

	int width, height, depth;
	string path;
	int start_frame, end_frame;



	width = depth = 100;
	height = 350; //(int)(width * 1.5f);
	path = GRIDY_RESULT_PATH;
	//start_frame = 1;
	end_frame = 500;


	if (argc >= 5)
	{
		width = atoi(argv[1]);
		height = atoi(argv[2]);
		depth = atoi(argv[3]);
		path = argv[4];

		if (argv[5]) end_frame = atoi(argv[5]);
	}


	GRIDY_SMOKE_3D mSMOKE;



	mSMOKE.Init(width, height, depth, path);
	mSMOKE.Simulate(1, end_frame);




	//mSMOKE.Init("");

	return 0;
}