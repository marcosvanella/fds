#!/usr/bin/env python3
"""Insert write-only REFDUMP hooks into a scratch copy of FDS main.f90 (unpatched otherwise). Usage: instrument.py <main.f90>"""
import re, sys
p = sys.argv[1]
s = open(p).read()
def rep(old, new, count=1):
    global s
    assert s.count(old) >= 1, old
    s = s.replace(old, new, count)
rep("USE OPENMP_FDS\n", "USE OPENMP_FDS\nUSE REFDUMP\n")
# init after set-up: right before the main loop
rep("MAIN_LOOP: DO\n", "CALL REFDUMP_INIT\n\nMAIN_LOOP: DO\n")
# predictor: viscosity + mass finite differences
rep("""      CALL INSERT_ALL_PARTICLES(T,NM)
      IF (.NOT.SOLID_PHASE_ONLY .AND. .NOT.FREEZE_VELOCITY) CALL COMPUTE_VISCOSITY(NM,APPLY_TO_ESTIMATED_VARIABLES=.FALSE.)
      CALL MASS_FINITE_DIFFERENCES(NM)
""", """      CALL INSERT_ALL_PARTICLES(T,NM)
      CALL REFDUMP_BEGIN('VISC_P',NM,T,DT)
      IF (.NOT.SOLID_PHASE_ONLY .AND. .NOT.FREEZE_VELOCITY) CALL COMPUTE_VISCOSITY(NM,APPLY_TO_ESTIMATED_VARIABLES=.FALSE.)
      CALL REFDUMP_END('VISC_P',NM)
      CALL REFDUMP_BEGIN('DENS_P',NM,T,DT)
      CALL MASS_FINITE_DIFFERENCES(NM)
""")
rep("""         CALL DENSITY(T,DT,NM)
         IF (LEVEL_SET_MODE>0) CALL LEVEL_SET_FIRESPREAD(T,DT,NM)
      ENDDO COMPUTE_DENSITY_LOOP""", """         CALL DENSITY(T,DT,NM)
         IF (FIRST_PASS) CALL REFDUMP_END('DENS_P',NM)
         IF (LEVEL_SET_MODE>0) CALL LEVEL_SET_FIRESPREAD(T,DT,NM)
      ENDDO COMPUTE_DENSITY_LOOP""")
rep("""            IF (.NOT.CYLINDRICAL) CALL VELOCITY_FLUX(T,DT,NM,APPLY_TO_ESTIMATED_VARIABLES=.FALSE.)
""", """            CALL REFDUMP_BEGIN('VFLUX_P',NM,T,DT)
            IF (.NOT.CYLINDRICAL) CALL VELOCITY_FLUX(T,DT,NM,APPLY_TO_ESTIMATED_VARIABLES=.FALSE.)
            CALL REFDUMP_END('VFLUX_P',NM)
""")
rep("""         CALL WALL_BC(T,DT,NM)
         IF (PARTICLE_DRAG) CALL PARTICLE_MOMENTUM_TRANSFER(NM)
         CALL DIVERGENCE_PART_1(T,DT,NM)
""", """         CALL WALL_BC(T,DT,NM)
         IF (PARTICLE_DRAG) CALL PARTICLE_MOMENTUM_TRANSFER(NM)
         CALL REFDUMP_BEGIN('DIV1_P',NM,T,DT)
         CALL DIVERGENCE_PART_1(T,DT,NM)
         CALL REFDUMP_END('DIV1_P',NM)
""")
rep("""         CALL DIVERGENCE_PART_2(DT,NM)
      ENDDO FINISH_DIVERGENCE_LOOP""", """         CALL REFDUMP_BEGIN('DIV2_P',NM,T,DT)
         CALL DIVERGENCE_PART_2(DT,NM)
         CALL REFDUMP_END('DIV2_P',NM)
      ENDDO FINISH_DIVERGENCE_LOOP""")
rep("""         CALL VELOCITY_PREDICTOR(T+DT,DT,DT_NEW,NM)
""", """         CALL REFDUMP_BEGIN('VPRED',NM,T+DT,DT)
         CALL VELOCITY_PREDICTOR(T+DT,DT,DT_NEW,NM)
         CALL REFDUMP_END('VPRED',NM,DT_NEW(NM),CHANGE_TIME_STEP_INDEX(NM))
""")
# corrector
rep("""      IF (.NOT.SOLID_PHASE_ONLY .AND. .NOT.FREEZE_VELOCITY) CALL COMPUTE_VISCOSITY(NM,APPLY_TO_ESTIMATED_VARIABLES=.TRUE.)
      CALL MASS_FINITE_DIFFERENCES(NM)
      CALL DENSITY(T,DT,NM)
""", """      CALL REFDUMP_BEGIN('VISC_C',NM,T,DT)
      IF (.NOT.SOLID_PHASE_ONLY .AND. .NOT.FREEZE_VELOCITY) CALL COMPUTE_VISCOSITY(NM,APPLY_TO_ESTIMATED_VARIABLES=.TRUE.)
      CALL REFDUMP_END('VISC_C',NM)
      CALL REFDUMP_BEGIN('DENS_C',NM,T,DT)
      CALL MASS_FINITE_DIFFERENCES(NM)
      CALL DENSITY(T,DT,NM)
      CALL REFDUMP_END('DENS_C',NM)
""")
rep("""         IF (.NOT.CYLINDRICAL) CALL VELOCITY_FLUX(T,DT,NM,APPLY_TO_ESTIMATED_VARIABLES=.TRUE.)
""", """         CALL REFDUMP_BEGIN('VFLUX_C',NM,T,DT)
         IF (.NOT.CYLINDRICAL) CALL VELOCITY_FLUX(T,DT,NM,APPLY_TO_ESTIMATED_VARIABLES=.TRUE.)
         CALL REFDUMP_END('VFLUX_C',NM)
""")
rep("""      CALL DIVERGENCE_PART_1(T,DT,NM)
""", """      CALL REFDUMP_BEGIN('DIV1_C',NM,T,DT)
      CALL DIVERGENCE_PART_1(T,DT,NM)
      CALL REFDUMP_END('DIV1_C',NM)
""", 2) if False else None
# second DIVERGENCE_PART_1 (corrector) is the one after the "COMPUTE_DIVERGENCE_2" area: find the last occurrence
i = s.rfind("      CALL DIVERGENCE_PART_1(T,DT,NM)\n")
assert i > 0
s = s[:i] + "      CALL REFDUMP_BEGIN('DIV1_C',NM,T,DT)\n      CALL DIVERGENCE_PART_1(T,DT,NM)\n      CALL REFDUMP_END('DIV1_C',NM)\n" + s[i+len("      CALL DIVERGENCE_PART_1(T,DT,NM)\n"):]
i = s.rfind("      CALL DIVERGENCE_PART_2(DT,NM)\n")
s = s[:i] + "      CALL REFDUMP_BEGIN('DIV2_C',NM,T,DT)\n      CALL DIVERGENCE_PART_2(DT,NM)\n      CALL REFDUMP_END('DIV2_C',NM)\n" + s[i+len("      CALL DIVERGENCE_PART_2(DT,NM)\n"):]
i = s.rfind("      CALL VELOCITY_CORRECTOR(T,DT,NM)\n")
s = s[:i] + "      CALL REFDUMP_BEGIN('VCORR',NM,T,DT)\n      CALL VELOCITY_CORRECTOR(T,DT,NM)\n      CALL REFDUMP_END('VCORR',NM)\n" + s[i+len("      CALL VELOCITY_CORRECTOR(T,DT,NM)\n"):]
# per-step record: T and DT used in this step (T before advance), DT_NEW(1), index; written right after "T = T + DT" in the corrector start
rep("   T = T + DT\n\n   ! Zero out energy and mass balance arrays", "   CALL REFDUMP_STEP(T,DT,DT_NEW(1),CHANGE_TIME_STEP_INDEX(1),0)\n   T = T + DT\n\n   ! Zero out energy and mass balance arrays")
rep("      IF (ANY(CHANGE_TIME_STEP_INDEX==-1)) THEN  ! If the time step was reduced, CYCLE CHANGE_TIME_STEP_LOOP",
    "      CALL REFDUMP_PASS(T,DT,DT_NEW(1),CHANGE_TIME_STEP_INDEX(1))\n      IF (ANY(CHANGE_TIME_STEP_INDEX==-1)) THEN  ! If the time step was reduced, CYCLE CHANGE_TIME_STEP_LOOP")
# finish before END_FDS
rep("ENDDO MAIN_LOOP\n", "ENDDO MAIN_LOOP\nCALL REFDUMP_FINISH\n")
open(p, 'w').write(s)
print("instrumented", p)
