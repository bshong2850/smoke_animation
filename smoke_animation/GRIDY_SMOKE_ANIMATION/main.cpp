#include <iostream>
#include "GRIDY_SMOKE_ANIMATION.h"

#include "../LIBRARY/GRIDY_CONFIG.h"

int main(int argc, char **argv)
{

	int width, height, depth;
	string path;
	string raymarching_path;

	int start_ani_frame, end_ani_frame, end_total_frame;


	width = 300;
	height = 220;
	depth = 90;
	path = GRIDY_RESULT_PATH;
	//raymarching_path = GRIDY_RAYMARCHING_PATH;
	start_ani_frame = 1;
	end_ani_frame = 16;
	end_total_frame = 16;


	if (argc >= 5) {
		width = atoi(argv[1]);
		height = atoi(argv[2]);
		depth = atoi(argv[3]);
		path = argv[4];

		if (argv[5]) end_ani_frame = atoi(argv[5]);
	}


	GRIDY_SMOKE_ANIMATED_MESH mSMOKE;



	//mSMOKE.Simulate(width, height, depth, path, raymarching_path, start_ani_frame, end_ani_frame, end_total_frame);
	mSMOKE.Simulate(width, height, depth, path, start_ani_frame, end_ani_frame, end_total_frame);




	//mSMOKE.Init("");

	return 0;
}