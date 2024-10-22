#include "buchberger.h"

int main()
{
	LeadMonomialsEnv env("K4_vc.txt",
					 "K4_vc_x0.txt",
				   false,
				   true,
					 14,
					 "VertexCover",
					 4,
					 2,
					 0.1);
	
	env.reset();
	env.step(0);
	return 0;
}
