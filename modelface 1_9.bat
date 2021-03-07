@echo off
color 4f
cls
echo -           -  - - - -  - - - -    - - -  -      -     - - -  - - -
echo - -       - -  -     -  -      -   -      -      -     -      -    -
echo -  -     -  -  -     -  -       -  -      -      -     -      -   -
echo -   -   -   -  -     -  -       -  - - -  -      -     - - -  - -
echo -    - -    -  -     -  -       -  -      -      -     -      -   -
echo -     -     -  -     -  -      -   -      -      -     -      -    -
echo -           -  - - - -  - - - -    - - -  - - -  - - - - - -  -     - && echo.
echo This is an interface to modeller 10.0 version && echo.
echo Please cite modeller as N. Eswar, M. A. Marti-Renom, B. Webb, 
echo M. S. Madhusudhan, D. Eramian, M. Shen, U. Pieper, A. Sali. 
echo Comparative Protein Structure Modeling With MODELLER. 
echo Current Protocols in Bioinformatics,echo John Wiley and Sons, 
echo Supplement 15, 5.6.1-5.6.30, 2006 && echo.
echo for more information about this API send email to^:
echo nava20ir@gmail.com
echo before any use modeller should be installed
echo DEPARTMENT OF MEDICINAL CHEMISTRY, FACULTY OF PHARMACY, Shiraz, IRAN && echo.
timeout 20
cls
set snp=0
set multi=n
color 0f
reg query HKLM\SOFTWARE\modelface /v mod
if %ERRORLEVEL%==1 (
goto pat
)
set KEY=HKEY_LOCAL_MACHINE\Software\modelface
for /f "tokens=2,*" %%a in ('reg query %KEY% /v vermol ^| findstr vermol') do (
set vermol=%%b
)
for /f "tokens=2,*" %%a in ('reg query %KEY% /v subvermol ^| findstr subvermol') do (
set subvermol=%%b
)
for /f "tokens=2,*" %%a in ('reg query %KEY% /v mod ^| findstr mod') do (
set mod=%%b
)
goto sett
:pat
::echo ### MODELLER PARAMETERS ### && echo.
::set /p passw=Please Enter the password:
::if not %passw%==shirazpharmacy (
::echo the password is wrong
::echo please send an email to nava20ir@gmail.com for the password
::pause
::exit
)
set vermol=10.0
set /p vermol=What is the version of your modeller x.x Enter for Default=10.0 :
set subvermol=0
set /p subvermol=What is the sub version of 10 default=0:
echo Enter the installation path of modeller?
set mod=C:\Program Files\Modeller%vermol%
set /p mod=default path is %mod% :
echo %mod%
reg add HKLM\Software\modelface /v vermol /t REG_SZ /d "%vermol%" 
reg add HKLM\Software\modelface /v subvermol /t REG_SZ /d "%subvermol%" 
reg add HKLM\Software\modelface /v mod /t REG_SZ /d "%mod%"
:sett
if not exist "%mod%"\lib\x86_64-w64\*.exe (
echo Modeller executables were not found
pause
goto pat
)
set vers=%vermol%
echo %subvermol%
set MODINSTALL10v%subvermol%=%mod%
set path=%mod%\lib\x86_64-w64;%path%
cls
:opt0
color 0f
echo ######## JOB MENU ######## && echo.
echo Type 1  and Enter To  Start Modeling with FASTA generation of template PDBs
echo Type 2  and Enter To  Run  alignment  section using Modeller %vermol% 
echo Type 3  and Enter For Model building  section using Modeller %vermol% 
echo Type 4  and Enter For Loop Refinement section using Modeller %vermol% 
echo Type 5  and Enter For Minimization  in vaccue using Modeller %vermol% 
echo Type 6  and Enter To  Get Energy and clashes  using Modeller %vermol% 
echo Type 7  and Enter For Modeling  Helix segment using Modeller %vermol% 
echo Type 8  and Enter For Repairing missing atoms using Modeller %vermol%
echo Type 9  and Enter To  Make Single point mutation in a PDB structure 
echo Type D  and Enter To  Delete this API setting parameters from registry
echo Type E  and Enter To  Exit && echo.
set task=1
set /p task=select one of the above operations:
goto opt%task%
:opt9
set snp=0
set /p snp=Please specify the number at which you would like to make mutant:
if %snp%==0  (
cls
echo No mutation point was defined && echo.
goto opt0
)
set /p multi=Is the native protien a mutli chain protein y/n:
set chk=A
:opt1
cls
if not exist *.pdb (
goto error
)
echo Generating FASTA files && echo.
FOR  %%i IN (*.pdb) DO ( 
echo %%~ni.fasta && echo.
del main%%~ni.txt
del %%~ni.fasta

findstr /B "ATOM" %%i|findstr /c:"CA" > main%%~ni.txt
setlocal enabledelayedexpansion
set fin=
set count=1
 FOR /F "tokens=4,5" %%S  IN (main%%~ni.txt) DO (
      set unk=0
      if %multi%==y (
       set ch=%%T
       set first=!ch:~0,1!
          if not !first! ==!chk! (
          set chk=!first!
          set fin=!fin!/
          )
      )   
  

if not !count!==%snp% (
       if %%S==ALA (
       set fin=!fin!A
       set unk=1
       )
       if %%S==ARG (
       set fin=!fin!R
       set unk=1
       ) 
       if %%S==ASN (
       set fin=!fin!N
       set unk=1
       )
       if %%S==ASP (
       set fin=!fin!D
       set unk=1
       )
       if %%S==GLU (
       set fin=!fin!E
       set unk=1
       )
       if %%S==GLN (
       set fin=!fin!Q
       set unk=1
       )
       if %%S==GLY (
       set fin=!fin!G
       set unk=1
       )
       if %%S==HIS (
       set fin=!fin!H
       set unk=1
       )
       if %%S==ILE (
       set fin=!fin!I
       set unk=1
       )
       if %%S==LEU (
       set fin=!fin!L
       set unk=1
       )
       if %%S==LYS (
       set fin=!fin!K
       set unk=1
       )
       if %%S==MET (
       set fin=!fin!M
       set unk=1
       )
       if %%S==PHE (
       set fin=!fin!F
       set unk=1
       )
       if %%S==PRO (
       set fin=!fin!P
       set unk=1
       )
       if %%S==SER (
       set fin=!fin!S
       set unk=1
       )
       if %%S==THR (
       set fin=!fin!T
       set unk=1
       )
       if %%S==TRP (
       set fin=!fin!W
       set unk=1
       )
       if %%S==TYR (
       set fin=!fin!Y
       set unk=1
       )
       if %%S==VAL (
       set fin=!fin!V
       set unk=1
       )
       if %%S==CYS (
       set fin=!fin!C
       set unk=1
       )
      
       if  !unk!==0 (
       set fin=!fin!.
       )


      
   ) else (
   echo.
   echo the residue at position !count! in native protein is %%S && echo.
   echo ALA=A  ARG=R   ASN=N  CYS=C
   echo ASP=D  GLU=E   GLN=Q  VAL=V   
   echo GLY=G  HIS=H   ILE=I  TYR=Y
   echo LEU=L  LYS=K   MET=M  TRP=W   
   echo PHE=F  PRO=P   SER=S  THR=T
   set /p mutant=Enter the one letter code of the residue at this site for mutant protein:
   set fin=!fin!!mutant!
   )

 set /a count=!count!+1
 )
echo !fin!
 if %snp%==0 (
 echo !fin! > %%~ni.fasta
 ) else (
 echo !fin! > seq.txt
 set snp=0
 set multi=n
 goto opt1
 )
del main%%~ni.txt
)
set proce=0
set /p proce=Press ENTER for alignment section or press 1 and Enter to exit:
if %proce%==1 goto opt0
:opt2
cls
color 3f
echo ### CLUSTAL SECTION ### && echo.
echo the files needed  are: && echo.
echo 1)A text file containing the FASTA of the unknown protein && echo.
echo 2)The fasta files of the template pdbs
pause
if not exist *.fasta (
cls
echo no fasta was found && echo.
color 0f
goto opt0
)
set k=seq.txt
:seq
dir/w *.txt
set /p k= Enter the unknown protein sequence file? Enter for Default=seq.txt:
if not exist %k% (
echo.
echo Error: no %k% does exist please try again
goto seq
)
echo %k%>seq_name.txt
set simmat=blosum62
set /p simmat=Enter the matrix name, blosum62 or as1, default is blosum62:
del structure.ali
del pdb_names.txt
cls
FOR  %%i IN (*.fasta) DO ( 
echo %%~ni>>pdb_names.txt
echo ^>P1;%%~ni>>structure.ali
echo structureX:%%~ni:first    :@ :end   : ::: :>>structure.ali
type %%~ni.fasta>>structure.ali
echo *>>structure.ali
)
echo ^>P1;%k% >temp.ali
echo sequence:%k%:first    :@ :end   : ::: :>>temp.ali
type %k%>>temp.ali
echo *>>temp.ali
echo from modeller import *>clustal.py
echo log.verbose()>>clustal.py
echo env = environ()>>clustal.py
echo env.libs.topology.read(file='$(LIB)/top_heav.lib')>>clustal.py
echo aln = alignment(env)>>clustal.py
echo aln.append(file='structure.ali', align_codes='all')>>clustal.py
echo aln_block = len(aln)>>clustal.py
echo aln.append(file='temp.ali', align_codes='%k%')>>clustal.py
echo aln.salign(rr_file='$(LIB)/%simmat%.sim.mat',  # Substitution matrix used>>clustal.py
echo            output='',>>clustal.py
echo            max_gap_length=20,>>clustal.py
echo            gap_function=True,              # If False then align2d not done>>clustal.py
echo            feature_weights=(1., 0., 0., 0., 0., 0.),>>clustal.py
echo            gap_penalties_1d=(-100, 0),>>clustal.py
echo            gap_penalties_2d=(3.5, 3.5, 3.5, 0.2, 4.0, 6.5, 2.0, 0.0, 0.0),>>clustal.py
echo            # d.p. score matrix>>clustal.py
echo            #write_weights=True, output_weights_file='salign.mtx'>>clustal.py
echo            similarity_flag=True)   # Ensuring that the dynamic programming>>clustal.py
echo                                    # matrix is not scaled to a difference matrix>>clustal.py
echo aln.write(file='alignment.ali', alignment_format='PIR')>>clustal.py
echo clustal.py has been generated
echo Ready for runing clustal
pause
echo please wait
mod%vers% clustal.py
del structure.ali
del temp.ali
if not exist alignment.ali (
type clustal.log
echo Error: no alignment.ali was generated Try again
goto opt2
)
set proce=0
set /p proce=Press ENTER for Model Building or press 1 and Enter to exit:
if %proce%==1 goto opt0
:opt3
cls
echo ### HOMOLOGY MODELLING SECTION ### && echo.
color 4f
echo the files needed  are: && echo.
echo 1)PDB(s) for the template proteins
echo 2)An alignment file with *.ali extension
echo 3)pdb_names.txt A text file containing the names of the template proteins
echo 4)seq.txt A text file containing the sequence fasta of unknown protein && echo.
pause
if not exist pdb_names.txt (
cls
echo pdb_names.txt not found && echo.
color 0f
goto opt0
)
set /p k=<seq_name.txt
set /p firstline=<pdb_names.txt
setlocal enabledelayedexpansion
set firstline='%firstline%'
FOR /F "skip=1" %%G IN (pdb_names.txt) DO (
set firstline=!firstline!,'%%G'
)
echo The pdb(s) %firstline% are being used as template(s) && echo.
set comman=automodel
set ifh=1
set num=5
set symetr=n
set disu=n
set md=n
set sss=
set /p ifh=Add hydrogens to the model? 1=No, 2=yes default=1:
if %ifh%==1 set comman=automodel
if %ifh%==2 set comman=allhmodel 
set /p num=How many models to build for the protein?Enter for Default=5:
set align=alignment.ali
:ali
@ echo on
dir/w *.ali
@echo off
set /p align=Enter the name of your aligmnent file?Enter for Default=alignment.ali:
if not exist %align% (
echo Error: no %align% does exist please try again
goto ali
)
echo when you have more than one chain, you can ask for the symmetry in the models
set /p symetr=Do you ask for symmetry in the models? y or n, default =n :
if %symetr%==y set /p ch=how many chains are there in your model?:
set /p disu=Is there any disulfide bridge? y or n:
set /p md=Perform MD minimization? y or n:
echo from modeller.automodel import *>model.py
if %symetr%==y echo def special_restraints^(self, aln^):>>model.py
if %symetr%==y (
FOR /L %%H IN (1,1,%ch%) DO (
set /p chname=Enter the %%H chain name 
set sss=!sss!s%%H,
echo  s%%H = selection^(self.chains['!chname!']^)>>model.py
)
)
if %symetr%==y echo  self.restraints.symmetry.append^(symmetry^(%sss% 1.0^)^)>>model.py
if %symetr%==y echo def user_after_single_model^(self^):>>model.py
if %symetr%==y echo  self.restraints.symmetry.report(1.0)>>model.py
if %disu%==y echo def special_patches(self, aln):>>model.py
if %disu%==y echo        self.patch(residue_type='DISU',>>model.py 
if %disu%==y set /p cysone=Enter the first  cys number:
if %disu%==y set /p cystwo=Enter the second cys number:
if %disu%==y echo                   residues=(self.residues['%cysone%'],self.residues['%cystwo%']))>>model.py
echo log.verbose()>>model.py
echo env = environ()>>model.py
echo env.io.atom_files_directory = './:../atom_files'>>model.py
echo #env.io.hetatm = True>>model.py                                                  
echo #env.io.water  = True>>model.py                                                  
echo a = %comman%(env, alnfile = '%align%', knowns = (%firstline%),>>model.py
echo               sequence = '%k%', assess_methods=(assess.DOPE, assess.GA341))>>model.py
echo a.starting_model= 1 >>model.py
echo a.ending_model  = %num% >>model.py
echo #a.deviation = 4.0>>model.py
if %md%==y echo a.md_level=refine.slow>>model.py
echo a.make()>>model.py
echo # Get a list of all successfully built models from a.outputs>>model.py
echo ok_models = filter(lambda x: x['failure'] is None, a.outputs)>>model.py
echo # Rank the models by DOPE score>>model.py
echo key = 'DOPE score'>>model.py
echo ok_models.sort(lambda a,b: cmp(a[key], b[key]))>>model.py
echo # Get top model>>model.py
echo m = ok_models[0]>>model.py
echo print "The best model : %%s (DOPE score %%.3f)" %% (m['name'], m[key])>>model.py
echo model.py has been generated
echo Ready for running MODELLER
pause
echo wait for the models to be generated && echo.
mod%vers% model.py
findstr /C:"The best model :" model.log >bestmodel.log
FOR /F "tokens=5" %%Q IN (bestmodel.log) DO set fff=%%Q
del bestmodel.log
FIND "The best model" model.log
set proce=0
set /p proce=Press ENTER for loop refinement or press 1 and Enter to exit:
if %proce%==1 goto opt0
:opt4
cls
echo ### LOOP REFINEMENT ### && echo.
color 2f
echo the files needed  are: && echo.
echo 1)The protein pdb file && echo.
pause
if not exist *.pdb (
goto error
)
echo from modeller import *>loop.py
echo from modeller.automodel import *>>loop.py
echo log.verbose()>>loop.py
echo env = environ()>>loop.py
echo env.io.atom_files_directory = './:../atom_files'>>loop.py
echo class MyLoop(loopmodel):>>loop.py
echo    def select_loop_atoms(self):>>loop.py
setlocal enabledelayedexpansion
set number=1
set /p number=How many loop regions are in the protien default=1?:
set /p looos=Enter the first loop START residue i.e 1:A :
set /p loooe=Enter the first loop END residue i.e 2:A :
echo           return selection(self.residue_range^('%looos%', '%loooe%'^)>>loop.py
if %number% GTR 1 (
FOR /L %%R IN (2,1,%number%) DO (
set /p looos=Enter the %%R th loop START residue i.e 1:A :
set /p loooe=Enter the %%R th loop END residue i.e 2:A :
echo ,self.residue_range^('!looos!', '!loooe!'^)>>loop.py
)
)
echo ^)>>loop.py
echo m = MyLoop(env,>>loop.py
set modeloop=%fff%
set /p modeloop=Enter the protein pdb file? default is %fff% :
echo          inimodel='%modeloop%', # initial model of the target>>loop.py
echo           sequence='looprefined', loop_assess_methods=assess.DOPE)               # code of the target>>loop.py
echo m.loop.starting_model= 1           # index of the first loop model>>loop.py
set numb=3
set /p numb=How many loop models to buid?press ENTER for default=3 :
echo m.loop.ending_model  = %numb%         # index of the last loop model>>loop.py
echo m.loop.md_level = refine.slow  # loop refinement method>>loop.py
echo m.make()>>loop.py
echo key = 'DOPE score'>>loop.py
echo loop.py has been generated
pause
echo wait
mod%vers% loop.py
echo                                      molpdf        DOPE
echo ------------------------------------------------------------
findstr /Birc:"looprefined.*.pdb" loop.log
pause
cls
color 0f
goto opt0
:opt5
cls
color 1f
echo ### Simulation ###
echo the files needed  are: && echo.
echo 1)The protein pdb file && echo.
pause
if not exist *.pdb (
goto error
)
:code
@ echo on
dir/w *.pdb
@echo off
set /p code=Enter the protein full name :
if not exist %code% (
echo.
echo Error: no %code% does exist please try again
goto code
)
set max_iter=20
set /p max_iter=Enter max iterations for first minimization step default=20:
set temp=300
set /p temp=Enter tempreature for MD step default=300:
set max_itermd=50
set /p max_itermd=Enter max iterations for MD step default=50:
set max_iter_sec=20
set /p max_iter_sec=Enter max iterations for second minimization step default=20:
set sav=10
set /p sav=Save protein coordinates every default=10:
set met=1
echo Choose the method press 1 for conjugate_gradients
echo                   press 2 for quasi_newton
set /p met=Enter the method:
echo %met%
echo from modeller import *>optimize.py
echo from modeller.scripts import complete_pdb>>optimize.py
echo from modeller.optimizers import conjugate_gradients, molecular_dynamics, quasi_newton, actions>>optimize.py
echo env = environ()>>optimize.py
echo env.io.atom_files_directory = ['../atom_files']>>optimize.py
echo env.edat.dynamic_sphere = True>>optimize.py
echo env.libs.topology.read(file='$(LIB)/top_heav.lib')>>optimize.py
echo env.libs.parameters.read(file='$(LIB)/par.lib')>>optimize.py
echo code = '%code%'>>optimize.py
echo mdl = complete_pdb(env, code)>>optimize.py
echo mdl.write(file=code+'.ini')>>optimize.py
echo atmsel = selection(mdl)>>optimize.py
echo mdl.restraints.make(atmsel, restraint_type='stereo', spline_on_site=False)>>optimize.py
echo mdl.restraints.write(file=code+'.rsr')>>optimize.py
echo mpdf = atmsel.energy()>>optimize.py
set max_at=100.0
set min_at=0.01
set /p min_at=Enter the minimum atom shift default=0.01
if %met%==1 set dd=conjugate_gradients
if %met%==2 set dd=quasi_newton
echo %dd%
echo cg=%dd%^(output='NO_REPORT', min_atom_shift=%min_at%^)>>optimize.py
echo md = molecular_dynamics(output='REPORT')>>optimize.py
echo trcfil = open(code+'.D00000001', 'w')>>optimize.py
echo cg.optimize(atmsel, max_iterations=%max_iter%, actions=actions.trace(5, trcfil))>>optimize.py
echo md.optimize(atmsel, temperature=%temp%, max_iterations=%max_itermd%,>>optimize.py
echo             actions=[actions.write_structure(%sav%, code+'.D9999%%04d_simulated.pdb'),>>optimize.py
echo                      actions.trace(%sav%, trcfil)])>>optimize.py
echo cg.optimize(atmsel, max_iterations=%max_iter_sec%,>>optimize.py
echo             actions=[actions.trace(5, trcfil)])>>optimize.py
echo mpdf = atmsel.energy()>>optimize.py
echo mdl.write(file=code+'.B')>>optimize.py
mod%vers% optimize.py
pause
cls
goto opt0
:opt6
cls
color 3f
echo ### Energy calculation ###
echo the files needed  are: && echo.
echo 1)The protein pdb file && echo.
pause
if not exist *.pdb (
goto error
)
:fff
@ echo on
dir/w *.pdb
@echo off
set /p fff=Enter your Protein full name:
if not exist %fff% (
echo.
echo Error: no %fff% does exist please try again
goto fff
)
echo from modeller import *>energy.py
echo from modeller.scripts import complete_pdb>>energy.py
echo env = environ()>>energy.py
echo env.io.atom_files_directory = ['../atom_files']>>energy.py
echo env.libs.topology.read(file='$(LIB)/top_heav.lib')>>energy.py
echo env.libs.parameters.read(file='$(LIB)/par.lib')>>energy.py
echo mdl = complete_pdb(env, "%fff%")>>energy.py
echo atmsel = selection(mdl)>>energy.py
echo mdl.restraints.make(atmsel, restraint_type='stereo',spline_on_site=False)>>energy.py
echo (molpdf, terms) = atmsel.energy(edat=energy_data(dynamic_sphere=True))>>energy.py
echo print("Bond energy is %%.3f" %% terms[physical.bond])>>energy.py
mod%vers% energy.py
findstr /c:"serious non-bonded atom clash" energy.log
findstr /Birc:"Bond energy is" energy.log
pause
cls
goto opt0
:opt7
cls
color 3f
echo ### Helix generator ###
echo the files needed  are: && echo.
echo 1)The sequence text file && echo.
pause
if not exist *.txt (
cls
echo no text file was found && echo.
color 0f
goto opt0
)
set filename=seq.txt
set bish=100
set ravesh=1
set optim=conjugate_gradients
:filename
@ echo on
dir/w *.txt
@echo off
set /p filename=Enter the file name default is seq.txt:
if not exist %filename% (
echo.
echo Error: no %filename% does exist please try again
goto filename
)
set /p bish=Enter the maximum iterations defaut=100:
set /p sequen=<%filename%
echo from modeller import *>helix.py
echo from modeller.optimizers import conjugate_gradients>>helix.py
echo # Set up environment>>helix.py
echo e = environ()>>helix.py
echo e.libs.topology.read('${LIB}/top_heav.lib')>>helix.py
echo e.libs.parameters.read('${LIB}/par.lib')>>helix.py
echo # Build an extended chain model from primary sequence, and write it out>>helix.py
echo m = model(e)>>helix.py
echo m.build_sequence('%sequen%')>>helix.py
echo m.write(file='extended-chain.pdb')>>helix.py
echo # Make stereochemical restraints on all atoms>>helix.py
echo allatoms = selection(m)>>helix.py
echo m.restraints.make(allatoms, restraint_type='STEREO', spline_on_site=False)>>helix.py
echo # Constrain all residues to be alpha-helical>>helix.py
echo # (Could also use m.residue_range() rather than m.residues here.)>>helix.py
echo m.restraints.add(secondary_structure.alpha(m.residues))>>helix.py
echo # Get an optimized structure with CG, and write it out>>helix.py
echo cg = %optim%()>>helix.py
echo cg.optimize(allatoms, max_iterations=%bish%)>>helix.py
echo m.write(file='%filename%.pdb')>>helix.py
echo helix.py was generated
pause
mod%vers% helix.py
cls
goto opt0
:opt8
cls
color 1f
echo ### Repair missing atom types ###
echo the files needed  are: && echo.
echo 1)The protein pdb file && echo.
pause
if not exist *.pdb (
goto error
)
cls
:codep
@ echo on
dir/w *.pdb
@echo off
set /p codep=Enter the protein full name :
if not exist %codep% (
echo Error: no %codep% does exist please try again
goto codep
)
set renu=1
set renm=False
set /p renu=Renumber the residues from beginning press[1]=yes [2]=no default=1:
if %renu%==2 set renm=True
echo from modeller import *>missatom.py
echo from modeller.scripts import complete_pdb>>missatom.py
echo env = environ()>>missatom.py
echo env.io.atom_files_directory = ['../atom_files']>>missatom.py
echo env.libs.topology.read(file='$(LIB)/top_heav.lib')>>missatom.py
echo env.libs.parameters.read(file='$(LIB)/par.lib')>>missatom.py
echo mdl = complete_pdb(env, "%codep%", transfer_res_num=%renm%)>>missatom.py
echo mdl.write(file='%codep%_repaired.pdb', model_format='PDB')>>missatom.py
echo missatom.py is generated
pause
mod%vers% missatom.py
pause
cls
goto opt0
:optE
Exit
:error
cls
echo no pdbs were found && echo.
color 0f
goto opt0
:optD
reg query HKLM\SOFTWARE\modelface /v mod
if not %ERRORLEVEL%==1 (
reg delete HKLM\SOFTWARE\modelface /f
cls
goto pat
)
