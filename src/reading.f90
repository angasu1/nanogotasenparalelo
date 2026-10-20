Module input_values
  use types      
  use global_variables
  use strings      
  use phys_cons
  use math_subroutines
  implicit none

  integer(ik)::nargs,ierror
  character(len=200)::warn=''
  character(len=200) :: linepru
  character(len=200),dimension(80) :: args
  character::singlecar
  character(len=2)::doublecar
  character(len=200) :: indir
  save

Contains


  subroutine reading_data
  character(len=200) :: format_string
  character(len=20)::ctemp
  integer(ik)::itemp,j


    open(file=trim(outdir)//'data.inp',unit=25,Status='OLD',Action='read',IOSTAT=ierror)
    open(file=trim(outdir)//'in_conf.xyz',unit=30,Status='OLD',Action='read',IOSTAT=ierror)
    if (ierror.ne.0) then
       write(*,*)'inputfile name data.inp does not exist'       
       stop
    else
    endif

    call title_reading
    format_string = "(A,A,A,I1,A,I1)"
    write (indir,format_string) '../results/',trim(title)
    indir=adjustl(trim(indir))
    write(*,*) 'dir=',trim(indir)
    call system('mkdir -p  '//trim(indir))
    call system('cp ../results/data.inp '//trim(indir))
    call system('cp ../results/in_conf.xyz '//trim(indir))

    call par_reading
    allocate(att(n_at+nhe))
    call reading_coords





      mtot=0.0_rk
      do j=1,n_at
       call props_assignment(att(j)%nom,j)
       if (j.le.n_at) mtot=mtot+att(j)%ma
      enddo
      mtot=mtot/nmon
      do j=1,nhe
       call props_assignment(att(n_at+j)%nom,n_at+j)
      enddo


      Call results_file_opening




     return
   end subroutine reading_data

   subroutine title_reading
    integer(ik)::i,j,k
    character::aa

    i=0    
    do    
        read(25,'(A)',IOSTAT=ierror) linepru
        aa=trim(adjustl(linepru))
        
        
        if (i.gt.1000) ierror=1
        if (ierror.ne.0) stop 'error reading data.inp in title section'
        i=i+1
        call parse(linepru,' ',args,nargs)
    
        if (aa.eq.'#') cycle
        if (aa.eq.'&') exit

       
    
        if (trim(args(1)).eq.'title') then
               title=trim(args(2))
        elseif (trim(args(1)).eq.'fixr') then
                fixran=.true.
        elseif (trim(args(1)).eq.'st') then
                if (trim(args(2)).eq.'dmc') then
                        simtyp=1
                elseif (trim(args(2)).eq.'dmcs') then
                        simtyp=2
                elseif (trim(args(2)).eq.'dvr') then
                        simtyp=3
                elseif (trim(args(2)).eq.'adiag') then
                        simtyp=4
                endif
        elseif (trim(args(1)).eq.'mt') then
                if (trim(args(2)).eq.'lm') then !Linear monomer
                        moltyp=1
                elseif (trim(args(2)).eq.'ld') then !Linear dimer
                        moltyp=2
                        nmon=2
                        initd=1.9_rk*ar2bo
                elseif (trim(args(2)).eq.'st') then  !Simmetric top
                        moltyp=3
                elseif (trim(args(2)).eq.'as') then   !Asimmetric
                        moltyp=4
                endif
        elseif (trim(args(1)).eq.'mn') then
                molname=trim(args(2))
                call check_name
                allocate(coor(3,n_atxmol))
        elseif (trim(args(1)).eq.'pt') then
                call value(args(2),ptyp,ierror)
                if (ptyp.eq.2) then     
                open(file=trim(potfilesdir)//'polinomios'//trim(molname)//'.dat' &
                        &,unit=1000,IOSTAT=ierror)
                open(file=trim(potfilesdir)//'potlam'//trim(molname)//'.dat' &
                        &,unit=1100,IOSTAT=ierror)
                call read_pot_files
                endif                      
        else
                write(*,*) 'argument ',args(1),' not understood'
                stop
        endif
    
    enddo    
   return
  end subroutine title_reading




   subroutine par_reading

        if (simtyp.le.2) then
              call dmc_par
        elseif (simtyp.eq.3) then
              call dvr_par
        elseif (simtyp.eq.4) then
              call ad_par
                write(*,*) 'simulation type(st) not understood'
                stop
        endif


   return
   end subroutine par_reading

   subroutine dmc_par
    integer(ik)::i,j

    i=0       
    do
       read(25,'(A)',IOSTAT=ierror) linepru
       i=i+1
       if (i.gt.1000) ierror=1
       if (ierror.ne.0) stop 'error reading data.inp'

       call parse(linepru,' ',args,nargs)
       if (args(1)(1:1).eq.'#') cycle
       if (args(1)(1:1).eq.'&') exit

          if (trim(args(1)).eq."nw") then
             call value(args(2),nw,ierror)   
          elseif (trim(args(1)).eq."nhe") then
             call value(args(2),nhe,ierror)    
          elseif (trim(args(1)).eq."nmon") then
             call value(args(2),nmon,ierror)    
                  n_at=nmon*n_atxmol  
          elseif (trim(args(1)).eq."state") then
             call value(args(2),state,ierror)    
          elseif (trim(args(1)).eq."node") then
             call value(args(2),node,ierror)    
          elseif (trim(args(1)).eq."stps") then
             call value(args(2),nstps,ierror)    
          elseif (trim(args(1)).eq."dtau") then
             call value(args(2),dtau,ierror)    
          elseif (trim(args(1)).eq."nruns") then
             call value(args(2),nruns,ierror)    
          elseif (trim(args(1)).eq."is") then
                  is=.true.
          elseif (trim(args(1)).eq."hh") then
               hhon=.false.
          elseif (trim(args(1)).eq."mh") then
               hmon=.false.
          elseif (trim(args(1)).eq."mm") then
               mmon=.false.
          elseif (trim(args(1)).eq."frotmol") then
             call value(args(2),frotmol,ierror)    
          elseif (trim(args(1)).eq."fdc") then
             call value(args(2),frothe,ierror)    
          elseif (trim(args(1)).eq."potterms") then
                  do j=1,20             
                     call value(args(1+j),pottermsp(j),ierror)
                     if (trim(args(j+1)).eq.'%') then
                     pottermsp(j)=0     
                     exit               
                     endif              
                     !ntermspot=ntermspot+1
                     ntermspot = pottermsp(j)
                  enddo                 
                     allocate (potterms(ntermspot))
                     allocate (potweig(ntermspot))
                     do j=1,ntermspot   
                        potterms(j)=pottermsp(j)
                     enddo              
          elseif (trim(args(1)).eq."weights") then
                  do j=1,ntermspot      
                     call value(args(1+j),potweig(j),ierror)
                     if (ierror.ne.0) stop 'number of weights different than termpots'
                  enddo                 
          else
             write(*,*) 'suboption ',trim(args(2)),' of DMC no recognized!!'       
             stop
          endif

    enddo


   return
   end subroutine dmc_par

   subroutine dvr_par
    integer(ik)::i

    i=0       
    do
       read(25,'(A)',IOSTAT=ierror) linepru
       i=i+1
       if (i.gt.1000) ierror=1
       if (ierror.ne.0) stop 'error reading data.inp'

       call parse(linepru,' ',args,nargs)
       if (args(1)(1:1).eq.'#') cycle
       if (args(1)(1:1).eq.'&') exit

          if (trim(args(1)).eq."ens") then
            !call value(args(2),enstyp,ierror)   
          elseif (trim(args(1)).eq."T0") then
            !call value(args(2),T0,ierror)    
          else
             write(*,*) 'suboption ',trim(args(1)),' of DVR no recognized!!'       
             stop
          endif

    enddo


   return
   end subroutine dvr_par

   subroutine ad_par
    integer(ik)::i

    i=0       
    do
       read(25,'(A)',IOSTAT=ierror) linepru
       i=i+1
       if (i.gt.1000) ierror=1
       if (ierror.ne.0) stop 'error reading data.inp'

       call parse(linepru,' ',args,nargs)
       if (args(1)(1:1).eq.'#') cycle
       if (args(1)(1:1).eq.'&') exit

          if (trim(args(1)).eq."ens") then
            !call value(args(2),enstyp,ierror)   
          elseif (trim(args(1)).eq."T0") then
            !call value(args(2),T0,ierror)    
          else
             write(*,*) 'suboption ',trim(args(1)),' of adiag no recognized!!'       
             stop
          endif

    enddo


   return
   end subroutine ad_par

   subroutine check_name

           Select case (molname)        
           case('hcl')                  
            n_atxmol=2                  
           case('hbr')                  
            n_atxmol=2                  
           case('hf')                   
            n_atxmol=2                  
           case('hfcl')                 
            n_atxmol=2                  
           case('hcn')                  
            n_atxmol=3                  
           case default                 
            write(*,*) 'molname ',trim(molname),' not recognized'
            stop                        
           end Select  
   
     return 
   end subroutine check_name

   subroutine reading_coords
           integer(ik)::i,j,nn
           
           read(30,*) nn
           if (nn.ne.n_atxmol) stop 'n_at wrong in in_conf.xyz'
           read(30,*)

           do i=1,nn
           read(30,*) att(i)%nom,coor(:,i)
           enddo
           requil=rad(coor(:,1),coor(:,2))
            
           close(30)

           if (nmon.gt.1) then
                   do i=2,nmon
                       do j=1,nn       
                        att((i-1)*nn+j)%nom=att(j)%nom
                       enddo
                   enddo
           endif

           do i=1,nhe
           att(n_at+i)%nom='He'
           enddo


    return       
   end subroutine reading_coords


      subroutine results_file_opening
      character(len=100)::cnhe,estado,ncam,cicorrida,ctt,cdtau,nombre,cis,crfb

        select case(simtyp)
        case (1) 
          write(cnhe,*) nhe
          write(estado,*) state
          write(ncam,*) nw
          write(ctt,*) nstps/10000
          write(cdtau,*) int(dtau)
          cnhe=adjustl(cnhe)
          estado=adjustl(estado)
          ncam=adjustl(ncam)
          ctt=adjustl(ctt)
          cdtau=adjustl(cdtau)
          if (is) then
           cis='is'
          else
           cis='nis'
          endif
          cis=adjustl(cis)

         !File opening
         nombre=trim(molname)//trim(estado)//'st'//&
               &trim(ncam)//'w'//trim(cnhe)//'he'//trim(ctt)//'tt'&
               &//trim(cdtau)//'dt'//trim(cis)//'.dat'
           write(*,*) nombre
           write(*,*) trim(indir)
          open (unit=100,file=trim(indir)//'/en'//trim(nombre))
          open (unit=250,file=trim(indir)//'/dist'//trim(nombre))
          open (unit=300,file=trim(indir)//'/enm'//trim(nombre))
          open (unit=400,file=trim(indir)//'/idist'//trim(nombre))
          open (unit=500,file=trim(indir)//'/config.xyz')
          open (unit=600,file=trim(indir)//'/output.dat')
       end Select

      end subroutine results_file_opening
      
      subroutine read_pot_files 
            integer(ik)                :: k, j

            ! determinar el número de datos en filas
            nLinea = 0            
            do     
               read(1000,*,iostat=ierror) temp
               if (ierror /= 0) exit
               nLinea = nLinea + 1
            end do                 

            rewind(1000)   

           ! Colocar memoria para la matriz        
           allocate(plam(nLinea, nColumna))
           allocate(vlam(nLinea, nColumna))
                                        
           ! Leer la matriz                        
           do k = 1, nLinea       
              read(1000,1300) (plam(k,j),j = 1, nColumna)
              read(1100,1300) (vlam(k,j),j = 1, nColumna)
           end do ! i = 1, nLinea 
      
           rmaxpot = vlam(k-1,1)



1300 format(*(f28.16,2x))                 
      end subroutine read_pot_files 
 End Module input_values

