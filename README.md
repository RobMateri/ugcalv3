This readme is for the inportance of knowing what my project is. I am currently updating the ugcalv2ud.f from its Fortran code to C++. This file shall explain the information of importance with all the current files present.

README-FORTRAN.md -> The original read me file for the Fortran version of ugcal (hereby known as ugcal.f). Kept as a momento sort of thing

ugcalv2ud.f -> The original ugcal program (known as ugcal.f). It is written in Fortran and due to this and it relying on old cern librarys to function, it requires some updating.

Makefile -> The makefile for ugcal.f so that it can compile. Unfortunately, it will not compile unless the cernlib is present which it will not as it has been discontinued many moons ago. The cernlib was used as it had a minimizer present in it that was used in ugcal.f.

ideal_result.txt -> The ideal result from running the ugcal program after it is compiled. This is being used to aide in figuring out bugs with the code along as the final goal

ugcalv3_01.c -> The basic upgraded template of my work with updating ugcal.f to a C++ written language. This will most not change as much as other versions as it is being used as a template. Note that this does compile but does not do its job.

ugcalv3_02.c -> Also known as ugcal.root, it will be using more ROOT based methods within in an attempt to make it more efficient and so that it may be used with ROOT instead of as a stand alone. It currently has no ROOT methods being used with it and will be worked on later

ugcalv3_03.c -> Also known as ugcal.brent, it uses a Brent minimizer which is written into the code. It currently compiles but only returns 0s for all numerical energy values (Not ideal)
