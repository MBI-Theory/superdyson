exit 1

Last updated: Jan. 27, 2026
============

IMPORTANT:
IMPORTANT: PLEASE READ THIS ENTIRE FILE BEFORE USING SUPERDYSON, OR
IMPORTANT: ASKING ANY QUESTIONS ABOUT ITS INSTALLATION OR USE.
IMPORTANT:

Property evaluation for general CI expansions. All the quantum chemistry data
(orbitals and CI coefficients - we can do our own 1-electron integrals now) are
taken from GAMESS-US. Unfortunately, the GAMESS interface is not fully
automated - it is a sequence of scripts. Please note that (a simple-ish) GAMESS
source code modification is required to obtain reliable results (see below).

At the moment, we support:

- overlap integrals, including displaced geometries
- transition dipoles
- many other property operators (see braket_operator in sd_core.f90 for
  the complete list).
- transition density matrices
- Dyson orbitals and "exchange"/"cradle" corrections

Contents:
---------
 1. GAMESS-US source code modifications
 2. Installation of superdyson
 3. Examples and test cases
 4. Input description
 5. Using GAMESS-US for preparing superdyson inputs.
 6. Simplified input
 7. Caveats

1. GAMESS-US source code modifications

GAMESS-US is very fond of adjusting the phases of the molecular orbitals.
This happens in many places, and the adjustments are triggered by many
different scenarios. Unfortunately, this means that the CI expansions printed
in the output are sometimes (but not always) for the determinants constructed
from the orbitals <i>different</i> from those printed in the output or punched
in the .dat file.

This is very bad news for any code which needs to use GAMESS-US wavefunctions
to calculate properties. The only reliable solution I am aware of is to punch
the MOs right at the point where GAMESS starts the CI calculation; unfortunately,
this requires a (very simple) modification of the GAMESS-US source.

What to do:

1a. Copy "ps_savemos.src" from the superdyson home directory to the GAMESS-US
    source directory (e.g. ~/gamess/source).
1b. Edit "gamess.src" to insert a call to ps_savemos right before the non-CIS
    CI packages are called (CIS is handled elesewhere). The exact location will
    depend on the version of GAMESS; in 15JUL2024R2P1 it is:

  3294       IF(CITYP.EQ.CIS .OR. CITYP.EQ.SFCIS) THEN
  3295          CALL WFNCIS
  3296          RETURN
  3297       END IF
  3298 C
  3299       CALL PS_SAVEMOS('Orbitals at the entry to CI calculation')
  3300 C
  3301 C        CARRY OUT A DETERMINANT BASED NEO FULL CI COMPUTATION
  3302 C

1c. Edit "cisgrd.src" to insert a call to ps_savemos_cis right after the phases
    of the molecular orbitals have been adjusted, like this:

   506       CALL DAREAD(IDAF,IODA,X(LVEC),NBF3,15,0)
   507       CALL STFASE(X(LVEC),NBF,NBF,NQMT)
   508       CALL DAREAD(IDAF,IODA,X(IENG),NBF,17,0)
   509 C
   510       CALL PS_SAVEMOS_CIS(X(LVEC),
   511      &                    'Orbitals at the entry to CI calculation')
   512 C
   513    90 CONTINUE

1d. (optional) Increase the number of significant digits in the CI
    coefficient printout. Usually, GAMESS-US only prints 7 digits after
    the decimal for the CI coefficients; this will limit the accuracy
    of the results you can get from superdyson. The print statements you
    are looking for are in aldeci.src, algnci.src, and ormas1.src. They
    look something like this:

      IF(MASWRK) WRITE(IW,'(4A,F10.7)') CONA(1:NACT+2),'|',
     *                                  CONB(1:NACT+2),'|  ',CI(IPOS)

    Replace "F10.7" with "F17.14". The similar change in cisgrd.src needs
    to be made to two of the FORMAT statements, replacing the F12.8 format
    by F22.18.

   7070 FORMAT(7X,'ALPHA',3X,I5,12X,I5,10X,F12.8)
   7072 FORMAT(7X,'BETA ',3X,I5,12X,I5,10X,F12.8)

    Please note that this modification will change the main GAMESS
    output. If you have any tools or visualization programs parsing
    this output, they may fail to work with the modified version. If
    you are not the only user of GAMESS-US on your system, you may
    want to keep the modified version separate.

1e. Change GAMESS-US build scripts to include the additional file.
    In compall, you will need to add something like this:

  538 ./comp zmatrx
  539 ./comp ps_savemos   # <- This is the new line
  540 #
  541 # gamess-base-end-section
  542 #


    In lked, something like this should work:

  1359                           modules_dft.o mod_vvos.o modf77_aambs.o utils_strings.o modules_common.o messages.o \
  1360                           fcidump.o ps_savemos.o)    # <- We've added ps_savemos.o here
  1361 if ($GMS_MSUAUTO == true) then
  1362    set STANDARD_GAMESS_OBJ4=(ccsd3aacgreorder.o ccsd3aacgsum.o ccsd3aacgt1A00.o \

1f. Compile and link GAMESS-US as described in GAMESS documentation.

2. Installation of superdyson

You will need a working Fortran compiler, BLAS, and LAPACK. If you would like
to run superdyson on multiple cores, your compiler and libraries must support
OpenMP. The libraries must be configured for single-threaded execution (they
will be called from a parallel region).

2a. Edit the Makefile to choose your compiler and libraries. There are sample
    command lines for gfortran (recommended) and Intel Fortran (likely to produce
    slower OpenMP code); however, these will likely need some adjustments.

2b. Run "make". If everything goes right, at the end you should have an
    executable binary, "sd_v2.x" in the current directory. If things go
    wrong, you should seek help from your local computer support. Unless
    you have a running executable already, there is very little the author
    of the code can help you with.

3. Examples and test cases

Successful compilation does not guarantee the <i>correct</i> compilation
(this is especially true for the Intel compiler; dodgy BLAS and LAPACK
libraries have also been known to pop up in some Linux distributions).

The next step is to run the test set. It is found in the "test" subdirectory.
In order to run all tests, switch to the test subdirectory, and execute the
command:

./run-all.sh

Depending on the speed of your computer, the test set may take several hours
to complete. The test script will also verify some key output from the code.
If it reports "OK:", the test is likely completed successfully.

Some false negatives are possible for very small results. For example, the
test case "n2-lowdin_prop_nosym.inp" will be frequently reported as "NOTOK:".
The exact result in this case is zero; however due to small round-off errors
very small non-zero values may be computed, triggering test failure.  If this
happens, you will need to compare the output with the reference output included
with the test set. For example (assuming you have vimdiff installed):

vimdiff n2-lowdin_prop_nosym.out n2-lowdin_prop_nosym.out_ref

Finally, if you need to re-run some or all tests (e.g. after recompiling or
modifying the code), you will need to remove all .out files from the test
directory. Otherwise, the run-all.sh script will simply try to compare the
existing outputs to the reference without re-running the calculations.

4. Input description

This section describes the "V2" input for superdyson. There is also the
legacy "V1" input, implemented as separate codes (superdyson.f90 and
tr_1rdm.f90). Unless you know what and why you are doing, stay away from
these two.

Superdyson expects a Fortran NAMELIST "SD_V2" on the standard input. All
secondary input is controlled from this namelist. All output is to the
standard output. The following keywords can appear in the namelist. The
default value is given after the equals sign.

 TASK = 'braket'

  Task to perform. Can be one of:

  'braket' - Evaluate integrals of the form <bra|ket> and <bra|op|ket>.

             If the bra and ket wavefunctions have the same number of
             electrons, these become overlap and 1-electron property integrals.

             If the ket wavefunction has one electron less than the bra
             wavefunction, <bra|ket> is the Dyson orbital. In the <bra|op|ket>
             integral, the property operator is taken as a sum over all
             coordinates of the ket wavefunction. If op = 'dipole', <bra|op|ket>
             is the "craddle" orbital of [PRA 80, 063411 (2009)].

  '1-rdm'  - Calculates McWeeny's 1-particle reduced density matrix between the
             bra and the ket wavefunctions. The 1-rdm is represented by its
             singular value decomposition. When bra and ket coincide, the
             decomposition is identical to the natural orbitals and their
             occupation numbers.

  '1-srdm' - Similar to '1-rdm', but calculates spin-resolved density matrix,
             which is represented by four spin blocks (AA, AB, BA, and BB).
             The '1-rdm' result is equivalent to the sum of the "AA" and "BB"
             blocks of '1-srdm', but keep in mind that the singular vectors
             may well be different!

 BRAKET_OPERATOR = 'dipole'

   Definition of the operator for TASK='braket'. Can be one of:

   'none'     - No property matrix elements will be computed. Property
                evaluation can be expensive, especially if Lowdin rules
                must be used.

   'kinetic'  - Kinetic energy operator [0.5 * (d^2/dx^2 + d^2/dy^2 + d^2/dz^2)]

   'dipole'   - Coordinate expectation [(x,y,z)]. This is a 3-component vector.
                If the true dipole is needed, please do not forget to multiply
                by the electron charge!

   'velocity' - Velocity operator [(d/dx,d/dy,d/dz)]. This is a 3-component
                vector.

   'r-d/dr'   - Mixed coordinate expectation - coordinate derivative operator.
                [r_i (d/d r_j)]. The i index runs fast. This is a 3x3 matrix.
                Angular momentum operator is the anti-covariant part of this
                operator.

   '1/r'      - Coulomb nuclear attraction [|r-R|^-1]. The probe postion R
                must be specified as OP_CENTRE.

   'r'        - Linear attraction [|r-R|]. The probe position R must be
                specified as OP_CENTRE.

   'r**2'     - Harmonic attraction [|r-R|^2]. The probe position R must be
                specified as OP_CENTRE.

   'gauss'    - Gaussian attraction [Exp[-alpha*|r-R|^2]. The probe postion
                R as in OP_CENTRE. The Gaussian exponent is in OP_GAUSS.

   'gauss r'  - Gaussian-damped linear attraction [|r-R|*Exp[-alpha*|r-R|^2].
                The probe postion and exponent are in OP_CENTRE and OP_GAUSS.

   Many other operators are available in the integral package, but are not
   currently accessible through the BRACKET_OPERATOR input.
   See import_gamess.f90 for the complete list of operators implemented.

 OP_CENTRE = 0.0, 0.0, 0.0

   Operator position. See BRACKET_OPERATOR above.

 OP_GAUSS = 1.0

   Operator exponent. See BRACKET_OPERATOR above.

 COMMENT = ' '

   Arbitrary string. It will be included in the output, but is not otherwise
   used in any way.

 NMOBRA = 0

   Number of single-particle orbitals (not electrons!) which appear in the
   bra wavefunction. In the GAMESS CI terms, this is the sum of the number
   of frozen, core, and active MOs.

 NMOKET = 0

   Number of single-particle orbitals in the ket wavefunction.

 NDETBRA = 0

   Number of determinants in the bra wavefunction.

 NDETKET = 0

   Number of determinants in the ket wavefunction.

 EPS_CDET = 1e-5_rk

   Cut-off for the product of determinant amplitudes. The default cut-off
   is propably too loose for most quantitative calculations.

 DONT_CHECK_S2 = .false.

   Disable calculation of the norm and S^2 expectation values for the
   input wavefunctions. The algorithm we use for calculating <S^2> is
   rather naive (and therefore slow). For simple overlaps, <S^2> may
   dominate the overall calculation. If you are sure that all inputs
   are correct, this is a simple way of speeding up things a bit.

 FILE_CBRA = ' '

   Name of the file containing molecular geometry, basis set, and MO
   coefficient specification. Superdyson accepts a subset of GAMESS-US
   punch files. Spefically:

   1. The $DATA section must be in C1 symmetry. The basis set must be
      specified explicitly (list of the exponents/contractions, rather
      than names. 'L' shells are however not supported).
   2. No $ECP section is allowed
   3. Only one $VEC section may be present. At least NMOBRA orbitals
      must be present in the $VEC section.

   Section 5 below gives an example of how this file may be prepared.

 FILE_CKET = ' '

   Name of the MO file for the ket wavefunction (see FILE_CBRA).

   Some or all of the geometry, basis set, and the MO coefficients can
   differ from those used for the bra wavefunction. However, if the
   single-particle orbitals of the ket coincide with the orbitals used
   by the bra, the more numerically efficient Slater rules will be
   used.

 FILE_DETBRA = ' '

   Name of the file containing the determinant list. The file must
   contain at least NDETBRA determinants, with each determinant
   described by the amplitude (a real number), followed by NMOBRA
   integers. The integers give the occupation of the single-particle
   orbitals in this determinant (+1: alpha-spin electron; -1: beta-spin
   electron; 2: both alpha- and beta-spin electrons).

   All determinants must have the same number of electrons.

   The sign of the amplitude assumes the alpha/beta string ordering
   for the spin-orbitals in the determinant (all alpha-spin occupied
   spin-orbitals before all beta-spin occupied spin-orbitals). This
   is the same convention used internally in GAMESS.

   Section 5 below gives an example of how this file may be prepared.

 FILE_DETKET = ' '

   Similar to FILE_DETBRA, but for the ket wavefunction.

 EPS_INTEGRAL = 1e-10

   Tolerance used to locate non-zero blocks of the overlap and
   property matrices. Has no effect if USE_SYMMETRY=.false.

 USE_SYMMETRY = .true.

   When .true., superdyson will use block structure of the overlap
   and property matrices to speed up calculations. The effects can
   be quite dramatic, especially for calculations of the reduced
   density matric.

   I am not aware of any cases where USE_SYMMETRY=.T. slows the
   calculation down or produces an incorrect result, so there is
   no reason (other than debugging) to turn it off.

 DETAIL_TIMERS = .false.

   Produce a very detailed timing report from low-level routines.
   The overhead of DETAIL_TIMERS=.T. is severe.

 SIGN_BRA = 1.0

   The bra wavefunction is multiplied by the factor SIGN_BRA
   (which does not have to be a unity). This is useful for producing
   globally-consistent property surfaces without modifying the
   actual wavefunctions.

 SIGN_KET = 1.0

   As SIGN_BRA, but for the ket wavefunction.

5. Using GAMESS-US for preparing superdyson inputs.

 Several example input files, covering all major classes of calculations
 in superdyson, are available in the test subdirectory. Preparing these
 inputs with the modified version of GAMESS-US (see section 1 above) is
 simple, but requires several steps.

 Here, we will go through those steps using an example of a Dyson orbital
 calculation for a CO (X^1\Sigma^+) -> CO+ (X^2\Sigma^+) ionization. All
 GAMESS-US inputs and reference outputs are in the gamess-example/
 subdirectory. The reference outputs and data files are in gamess-example/ref_out

 5a. We begin with a closed-shell Hartree-Fock calculation for the neutral,
     to prepare starting orbitals for CASSCF calculations in the later
     steps. The input file is "co-rhf.inp". Some key input keywords here:

     ITOL=40 ICUT=20 <- Requests high-accuracy integrals. The defaults are not
                        sufficiently accurate in the asymptotic region.
     EXTFIL=.TRUE.   <- Requests that the basis should be read from an external
                        file (specified by the EXTBAS environment variable; see
                        gms-files.csh in the GAMESS-US installation directory).
                        Even if you want to use one of the basis sets GAMESS
                        knows about, you MUST nonetheless read it from an
                        external file to make sure that superdyson can parse
                        GAMESS output files. You can use http://bse.pnl.gov/
                        to generate the basis set files.
     GBASIS=APVDZ    <- Refers to a basis set in the "basis" file. This is
                        the aug-cc-pVDZ basis. It is likely too small for any
                        production work, but it makes everything fast ...

     Feel free to look up the rest of the keywords in the GAMESS-US manual.

     In addition to the GAMESS-US output, we will also need the PUNCH (.dat)
     file. It is usually found in ~/scr/ directory.

 5b. Next, we perform a CASSCF calculation for the CO(1+) cation. The input
     file is "cop-cas.inp". Again, some key parts of the input:

     GUESS=MOREAD    <- Fetch guess MOS from the $VEC section. We have
                        copied these orbitals from the .dat file produces
                        by the RHF calculation in step 5a above.
     NORB=46         <- Number of the MOs. Taken from the RHF output file
     CISTEP=GMCCI    <- We are using the GMCCI CI program, which allows
                        us to average over multiple irreps (so that we can
                        still benefit from the molecular symmetry in CASSCF).
                        It is however perfectly fine to use any CI program
                        in this step.
     KSTATE(1)=1,1,1,1,1,1,1
                     <- We are averaging our solutions over several state
                        of the cation (X^2\Sigma^+, A^2\Pi, B^2\Sigma^+,
                        D^2\Pi, and C^2\Sigma^+ in our case). Since we are
                        only interested in the X state here, this is not
                        really necessary.

    Again, in addition to the output we will need the punch (.dat) file from
    ~/scr/

 5c. Then, we perform a MR-CI calculation for the ground state of the CO(1+)
     cation. The input file is cop-cas-ci_a1.inp. Some key parts of the input:

     CITYP=ORMAS     <- We will use ORMAS CI program. You must choose either
                        ALDET, GENCI, or ORMAS here. It is likely that other
                        CI programs in GAMESS may also work, our conversion
                        script does not now how to handle the output.
     PURIFY=.FALSE.  <- We need to disable parts of GAMESS which modify input
                        MOs and may trigger a phase adjustment
     TOLE=0.0        <- ... and another possible orbital modification
     TOLZ=0.0        <- ... and another ...
     IPURFY=0        <- ... and another ...
     CUTTRF=1E-12    <- The default cut-off for direct integral transformation
                        is too loose when diffuse orbitals are present.
     PRTTOL=-1       <- We need a complete determonant list, not just the bits
                        with large amplitudes.

     $VEC            <- Please note that we have replaced the $VEC section with
                        the optimized CASSCF orbitals from step 5b above.
                        In principle, one case use any other set of single-particle
                        orbitals here, not necessarily the set coming from an
                        MCSCF calculation. This will affect the accuracy of the
                        results, but the proceduce for preparing superdyson inputs
                        remains the same.

    Although we are using the ORMAS program here, the active space we specify
    is actually the same CAS used in the CASSCF step. Therefore, the state
    energies must match the CASSCF output (check!).

    As usual, we need both the output file and the .dat file from ~/scr/
    Note that the main output (which now contains the complete determinant list)
    may be rather large.

 5d. We can now prepare the .mos and .dets files for the superdyson ket wavefunction.
     For the .mos file:

     5d1. Edit the $DATA section to switch to the C1 symmetry (see GAMESS manual
          for the $DATA section if you are not sure how).
     5d2. Make sure that the $VEC section immediately following the $DATA section
          has a commend line "ALPHA MOS: Orbitals at the entry to CI calculation"
     5d3. Delete everything after the $END matching the $VEC.
     5d4. You should now have something looking like ref_out/cop-cas-ci_a1.mos

     5d5. To construct the determinant list, you can use the command:

          ../det_conv.awk state=1 cop-cas-ci_a1.out > cop-cas-ci_a1.dets

          Parameter "state=1" requests the first state in the CI output. Please
          keep in mind that the determinant list may be rather large!

 5e. Next, we need to calculate the wavefunction for the neutral state. The steps
     are the same as in 5b, 5c, and 5d above. The inputs are:

      con-cas.inp        <- CASCF input
      con-cas-ci_a1.inp  <- MR-CI input

 5f. You now have everything you need to calculate the Dyson orbital. The
     example input file for superdyson is in "co_dyson.inp"

6. Simplified input

 In many cases, the full superdyson interface described above may be unneccessary
 (not to mention cumbersome and error-prone). A simplified interface is provided
 by the "overlap.sh" script. The calling sequence is:

   overlap.sh left_out left_state right_out right_state [task] [operator]

   {left,right}_out   should specify GAMESS .dat file containing the orbitals
                      and optionally the matching .out and/or .inp files
   {left,right}_state may be an integer state index, in which case GAMESS .out
                      specofoed by {left,right}_out must exist. Alternatively,
                      it may also give a determinant list in the internal
                      superdyson format.

   Optional task and operator are as defined in superdyson; see above

   All files may be compressed with bzip2.

 The script is able to handle inputs in the c1, d2, cs, c2v, c2h, and d2h point
 groups, provided that the default setting is used. It will recognize ALDET, GENCI,
 and ORMAS wavefunctions.

7. Caveats

 There is very little error and consistency checking in superdyson. It basically
 assumes that you know what you are doing. If you don't, the results may be
 disastrously wrong!

 Only the determinantal CIS wavefunctions are currently supported. Spin-flip CIS
 and CSF CIS are not.
