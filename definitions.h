#define  PHYSICS                        MHD
#define  DIMENSIONS                     1
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  HANCOCK
#define  NTRACER                        0
#define  PARTICLES                      PARTICLES_CR
#define  USER_DEF_PARAMETERS            4

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   NO
#define  BACKGROUND_FIELD               NO
#define  AMBIPOLAR_DIFFUSION            NO
#define  RESISTIVITY                    NO
#define  HALL_MHD                       NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  RADIATION                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  USER_VNOISE                    0
#define  USER_NCR_by_NGAS               1
#define  USER_UCR                       2
#define  USER_VBOOST                    3


/* [Beg] user-defined constants (do not change this line) */

#define  CT_EMF_AVERAGE                 ARITHMETIC
#define  LIMITER                        MC_LIM
#define  PARTICLES_CR_C                 50.
#define  PARTICLES_CR_E_MC              1.
#define  PARTICLES_CR_E_MC_GAS          1.
#define  PARTICLES_CR_NSUB              4
#define  PARTICLES_DEPOSIT              INTEGER
#define  PRIMITIVE_HANCOCK              FALSE
#define  VTK_TIME_INFO                  TRUE

/* [End] user-defined constants (do not change this line) */
