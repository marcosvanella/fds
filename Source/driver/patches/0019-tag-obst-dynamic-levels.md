# 0019 --tag-obst on dynamic levels (Source/regrid_transport/RegridAmrCore.cpp:156)
Status: **apply** (not committed; Role 2/3 file).
The --tag-obst block used the static-footprint guard `lev + 1 <= h_.top`; with level-0 &MESH only h_.top = 0, so the OBST tags were never set
and the pyramid stayed on level 0. Guard is now `lev < maxLevel()`.
