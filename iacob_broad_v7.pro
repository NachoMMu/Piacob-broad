@sub_compile.pro

pro iacob_broad,star,resol,d,line=line,cut=cut,noproc=noproc,resul=resul,ps=ps,$
	 xl=xl,xyc=xyc,xsnr=xsnr,vs_user=vs_user,vm_user=vm_user,$
	 fits=fits,vfts=vfts,nomads=nomads,$
	 gauss=gauss,rerec=rerec,vr=vr,em=em,beta=beta,sigma_FT0=sigma_FT0

; *******************************************************************************************************
;                                          VERY IMPORTANT !!! 
; *******************************************************************************************************
; (1) Indicate below the name of your ps-file viewer
; (2) Indicate below the path to the directory where iacob_broad is stored
; *******************************************************************************************************
 
  ps_viewer='okular' 
  
  home='/media/AURYN2/idlpro/VSINI/iacob_broad_public/'

; *******************************************************************************************************
;    Please, refer to Simon-Diaz & Herrero (2014, A&A 562, 135) in any publication using iacob-broad 
; *******************************************************************************************************

; =======================================================================================================
; 		        iacob_broad procedure (S. Simon-Diaz, v7, 2014/11/01)
; =======================================================================================================
;
; -------------------------------------------------------------------------------------------------------
; --- INPUT ---------------------------------------------------------------------------------------------
; -------------------------------------------------------------------------------------------------------
;
; star  	    : Short identification of the spectrum file (to be used in the output)
; resol 	    : Resolving power of the spectrum
;
; -------------------------------------------------------------------------------------------------------
; --- OUTPUT --------------------------------------------------------------------------------------------
; -------------------------------------------------------------------------------------------------------
;
; d		    : Structure with the summary of results
;
; In addition, the program creates 2 output xdr-files, which can be easily restored in IDL: 
;
;   e.g. HD37042_SiIII4552.xdr
;   e.g. HD37042_SiIII4552_resul.xdr or HD37042_SiIII4552_resul_vg.xdr (gauss=1)
;
; and a ps-file with a graphical summary of results
;
; See the manual for the description of what is stored in the xdr-files
;
; -------------------------------------------------------------------------------------------------------
; --- OPTIONAL KEYWORDS ---------------------------------------------------------------------------------
; -------------------------------------------------------------------------------------------------------
;  line=line	    : Line to be analyzed (selected from the lines indicated in iacob_broad_lines.dat)
;  cut=cut	    : The program exits after pre-processing of the line (SNR, spectra window, 
;		    : identification of the zero of the FT)
;  noproc=1	    : The observed line profile will not be pre-procesed (only valid if done before) 
;  resul=1	    : Pre-procesing of the line profile and GOF computation will not be done
;		    : (only valid if done before)
;  ps=1  	    : To open the .ps file
;
;  fits=1	    : To be used for IRAF-FITS spectra (see iacob_broad_leesp.pro)
;  vfts=1	    : To be used for VFTS-FITS spectra (see iacob_broad_leesp.pro) 
;  nomads=1	    : To be used for NoMaDS-FITS (or CAFE-BEANS-FITS) spectra (see iacob_broad_leesp.pro)
;
;  rerect=1	    : Renormalization option
;  vr=vr 	    : vrad correction option
;  em=1  	    : Emission line option
;
;  xl=[xa,xb]	    : Fixed spectral window used for the analisis of a give line profile 
;  xyc=[xc,yc]	    : Lambda of the core (xc) of the line and flux of the adjacent continuum (yc)
;  xsnr=[x1,x2]     : Fixed spectral window used for the SNR computation of the adjacent continuum
;
;  gauss=1	    : uses a Gaussian definition of the extra-broadening (default: RT)
;  vs_user, vm_user : Values of vsini and Theta_RT provide by the user for the final plots
;  beta=beta	    : limb darkening coeficient (default: 1.5 --> epsilon=0.6)
;  sigma_FT0	    : Frequecy of the first zero of the FT (default: 0.660)
; -------------------------------------------------------------------------------------------------------
;
; =======================================================================================================


; --------------------------------------------------------------------------------------------------------
;  Read the list of lines available (IF keyword 'line' is activated)
;  ==> New lines can be added to the iacob_broad_lines.txt file
; --------------------------------------------------------------------------------------------------------

  readcol,home+'iacob_broad_lines.dat',line_list,lam_list,FORMAT='(A,F)',/silent

; --------------------------------------------------------------------------------------------------------
; Initial window used for pre-processing of the spectrum line (IF keyword 'line' is activated)
; --------------------------------------------------------------------------------------------------------

  IF KEYWORD_SET(line) THEN BEGIN

    ll=where(line_list EQ line)

    lam0=lam_list(ll(0))    

    xsh=20.  &   xr0=lam0+xsh*[-1.,1]

    line_lab='_'+line
  
  END ELSE BEGIN

    line=''  &  line_lab=''
  
  END
  
; --------------------------------------------------------------------------------------------------------
; Minimum vsini and Theta_RT considered (FT plot and GOF part)
; --------------------------------------------------------------------------------------------------------

  vmin          = 0.25*2.99d5/resol  

; --------------------------------------------------------------------------------------------------------
; Maximum vsini and Theta_RT considered, relative to vsini(FT) (GOF part)
; --------------------------------------------------------------------------------------------------------
; ==> The GOF computation for vsini and vmac will run from vmin to fact_vs*vsini(zero FT) and 
;     fact_vrt*vsini(zero FT), respectively
;
; * NOTE: If a HeI is analyzed and an intrinsic delta-profile is assumed, a large fact_vrt value 
;         may be needed
; --------------------------------------------------------------------------------------------------------

  fact_vs       = 2.0   
  fact_vrt      = 2.5   

  IF STRPOS(strupcase(line),'HEI') GE 0 THEN fact_vrt=4.5 

; --------------------------------------------------------------------------------------------------------
; Step in vsini and Theta_RT for the GOF part (relative to vsini(FT), i.e., step_v*vsini(zero-FT))
; --------------------------------------------------------------------------------------------------------

  step_v        = 0.05  

; --------------------------------------------------------------------------------------------------------
; Range and step in EW (GOF part)
; --------------------------------------------------------------------------------------------------------

  ew_step       = 0.05  ;  Step in EW (relative to the EW provided by FT(sigma=0))

  n_ew_down     = 0.    ;  Number of steps with EW < EW0
  n_ew_up       = 0.    ;  Number of steps with EW > EW0

  n_ew_down_cut = 0.    ;  Number of steps with EW < EW0 (for the case when the line is clipped)
  n_ew_up_cut   = 3.    ;  Number of steps with EW > EW0 (for the case when the line is clipped)


; ========================================================================================================
; ========================================================================================================
;                  DO NOT MODIFY ANYTHING BELOW THIS POINT (OR DO IT AT YOUR OWN RISK!!!!)
; ========================================================================================================
; ========================================================================================================

  device,retain=2,true=24,decomposed=0
  quiet=1

  !p.multi = 0
  !p.thick = 1
  ilupal

  IF NOT KEYWORD_SET(vs_user)  THEN vs_user    = 0
  IF NOT KEYWORD_SET(vm_user)  THEN vm_user    = 0
  IF     KEYWORD_SET(gauss)    THEN labmac     = '_vg' ELSE labmac = ''
  IF     KEYWORD_SET(resul)    THEN noproc     = 1

  clip = 0

  ss=findfile(star+line_lab+'.xdr')
    
; ========================================================================================= BEGIN noproc=1
  IF NOT KEYWORD_SET(noproc) OR ss(0) EQ '' THEN BEGIN        ; < ========================= BEGIN noproc=1
; ========================================================================================= BEGIN noproc=1

    resul=0

    print,''
    IF ss(0) EQ '' THEN print,star+line_lab+'.xdr file not found ... lets create it ...'
  
  ; ------------------------------------------------------------------------------------------------------
  ; PREPARE THE LINE TO ANALYZE 
  ; ------------------------------------------------------------------------------------------------------

  ; ------------------------------------------------------------------------------------------------------
  ; --- Search and read the spectrum ---------------------------------------------------------------------
  ; ------------------------------------------------------------------------------------------------------

    iacob_broad_leesp,star,w,f,fits=fits,vfts=vfts,nomads=nomads,vr=vr

  ; ------------------------------------------------------------------------------------------------------
  ; --- Initial region where the line will be found ------------------------------------------------------
  ; ------------------------------------------------------------------------------------------------------

    IF line EQ '' AND NOT KEYWORD_SET(xl) THEN BEGIN

       plot,w,f,title='Select the region where the line is located',/ysty,/xsty
       plot_lin,lam_list
      
      xhair,x1,y & xhair,x2,y

      xr0=minmax([x1,x2])

    END
   
    tt=where(w GE xr0(0) AND w LE xr0(1))

    w0=w(tt) & f0=f(tt)

    yr0=minmax(f0)

    !y.range=[0.9*yr0(0),1.1*yr0(1) < 1.5]

  ; ------------------------------------------------------------------------------------------------------
  ; --- Local renormalization ----------------------------------------------------------------------------
  ; ------------------------------------------------------------------------------------------------------

    IF KEYWORD_SET(rerec) THEN BEGIN

      plot,w0,f0,title='Normalization: Select two continuum regions',/ysty,/xsty
      plot_lin,lam_list

      xhair,x1,y & xhair,x2,y & xr1=minmax([x1,x2])
      xhair,x3,y & xhair,x4,y & xr3=minmax([x3,x4]) 

      tt=where((w0 GE xr1(0) AND w0 LE xr1(1)) OR (w0 GE xr3(0) AND w0 LE xr3(1)))

      b=poly_fit(w0(tt),f0(tt),1)
      fr=poly(w0,b)

      f0=f0/fr

      yr0=minmax(f0)

      !y.range=[0.9*yr0(0),1.1*yr0(1) < 1.5]

    END

  ; ------------------------------------------------------------------------------------------------------
  ; --- Compute the SNR ----------------------------------------------------------------------------------
  ; ------------------------------------------------------------------------------------------------------

    IF NOT KEYWORD_SET(xsnr) THEN BEGIN
    
      plot,w0,f0,title='Mark a region to compute the SNR',/ysty,/xsty
      plot_lin,lam_list

      xhair,x1,y & xhair,x2,y

    END ELSE BEGIN
   
      x1=xsnr(0) & x2=xsnr(1)  
    
    END
    
    xr1=minmax([x1,x2])
    
    tt=where(w0 GE xr1(0) and w0 LE xr1(1))

    sigma=stddev(f0(tt)) & snr=1./sigma

    print,'x(SNR)= ',x1,x2
    print,'SNR = ',snr,FORMAT='(A6,I4)'
    print,''
 
  ; ------------------------------------------------------------------------------------------------------
  ; --- Indicate the limits of the line (for the FT and GOF computations) --------------------------------
  ; ------------------------------------------------------------------------------------------------------

    IF NOT KEYWORD_SET(xl) THEN BEGIN
    
      plot,w0,f0,title='Mark the limits of the line',/ysty,/xsty
      oplot,w0(tt),f0(tt),co=2
      plot_lin,lam_list

      xhair,x1,y & xhair,x2,y

    END ELSE BEGIN
  
      x1=xl(0) & x2=xl(1)
        
    END
    
    xr1=minmax([x1,x2])
    
    tt=where(w0 GE xr1(0) and w0 LE xr1(1))

    wx0=w0(tt) & fx0=f0(tt)

    !x.range=[min(wx0),max(wx0)]+0.1*(max(wx0)-min(wx0))*[-1,1]
    !y.range=[0.98*min(fx0),1.02*max(fx0) < 1.5]

  ; ------------------------------------------------------------------------------------------------------
  ; --- Indicate the continuum and core of the line ------------------------------------------------------
  ; ------------------------------------------------------------------------------------------------------

    IF NOT KEYWORD_SET(xyc) THEN BEGIN
    
      plot,w0,f0,/xsty,/ysty,title='Mark the continuum & the core of the line (just 1 click!)'
      oplot,wx0,fx0,co=2
      plot_lin,lam_list

      xhair,xc,ycont
 
    END ELSE BEGIN
    
      xc=xyc(0) & ycont=xyc(1)

    END

    fx0 = fx0 - ycont+1.
    f0  = f0  - ycont+1.

  PRINT,'xline =',x1,x2,xc,ycont
   
  ; ------------------------------------------------------------------------------------------------------
  ; --- Clip part of the line? ---------------------------------------------------------------------------
  ; ------------------------------------------------------------------------------------------------------

    IF NOT KEYWORD_SET(xl) THEN BEGIN
    
      a=1
 
      plot,w0,f0,/xsty,/ysty,title='Do you want to clip part of the line? (see terminal)'
      oplot,wx0,fx0,co=2,psym=2
      plot_lin,lam_list
 
      WHILE a DO BEGIN
 
        PRINT,'Do you want to clip part of the line? (0/1)'
        READ,a

        IF a EQ 1 THEN BEGIN

          plot,w0,f0,/xsty,/ysty,title='Mark the region of the line you want to clip'
          oplot,wx0,fx0,co=2,psym=2
          plot_lin,lam_list

          clip=1

          xhair,x1,y & xhair,x2,y

          xr1=minmax([x1,x2])
	  
          ; --- Linear interpolation (maybe consider a better option in the future?) ---
	
          tt_cut=where(wx0 GT xr1(0) and wx0 LT xr1(1))

          tt1=max(where(wx0 LE xr1(0))) & wx1=wx0(tt1) & fx1=fx0(tt1)
          tt2=min(where(wx0 GE xr1(1))) & wx2=wx0(tt2) & fx2=fx0(tt2)

          fx0(tt_cut)=(fx2-fx1)/(wx2-wx1)*(wx0(tt_cut)-wx1)+fx1

          plot,w0,f0,/xsty,/ysty,title='Mark the region of the line you want to clip'  
          oplot,wx0,fx0,co=2,psym=2
          plot_lin,lam_list

       END

     END

   END
   
   !x.range=0
   !y.range=0

  ; ------------------------------------------------------------------------------------------------------
  ; --- FT of the line -----------------------------------------------------------------------------------
  ; ------------------------------------------------------------------------------------------------------

    vfft,wx0,fx0,ycont,xc,base,factor,res,nr,max_sigma_obs,ew,nyq=1,sigma_FT0=sigma_FT0

    xfou_obs=1./(base*factor)
    yfou_obs=alog10(res(0:nr-1))
    
    a=0
    
    FOR i=0,1 DO BEGIN
    
      IF i EQ 0 THEN tit='Fourier transform, do you want to modify max(vsini) in the plot? (see terminal)'
      IF i EQ 1 THEN tit='Fourier transform, mark the position of the first zero'
    
      IF a EQ 0 THEN vft_max=450. ELSE vft_max=a
    
      plot,xfou_obs,yfou_obs,xtitle='v sini (km s!u-1!n)',$
           ytitle='ABS ( AMPLITUDE )',xr=[vft_max,0.5],/xstyle,$
           yr=[-4.,0.],/ystyle,title=tit
      oplot,[max_sigma_obs,max_sigma_obs],[-8,-0.7],line=4
      xyouts,max_sigma_obs,-0.6,'Nyq.',charsize=0.7,align=0,orientation=90

      IF i EQ 0 THEN BEGIN
        PRINT,'Do you want to modify max(vsini) in the FT-plot? Please, provide value (0: NO)'
        READ,a
      END
    
    END
    
    xhair,vft,y
    print,'vsini(FT) = ',vft,' km/s',FORMAT='(A12,F5.1,A5)'
    print,'sigma(FT) = ',1.d0/(factor*vft),' s-1',FORMAT='(A12,F8.3,A5)'

    ew=1.e3*ew

    PRINT,'EW (FT) = ',ew,' mA',FORMAT='(A10,F6.1,A3)'
    PRINT,''
       
    IF clip EQ 1 THEN BEGIN

      ; --------------------------------------------------------------------------------------------------
      ; Once the FT of the line is computed, the clipped part of the line is eliminated for 
      ; the GOF computation
      ; --------------------------------------------------------------------------------------------------
 
      fx0(tt_cut)=0

      tt=where(fx0 GT 0)
      wx0=wx0(tt)
      fx0=fx0(tt)
 
    END
 
    save,file=star+line_lab+'.xdr',w0,f0,wx0,fx0,xc,ew,snr,vft,xfou_obs,yfou_obs,max_sigma_obs,clip

    IF KEYWORD_SET(cut) THEN RETURN

; ========================================================================================== END noproc=1
  END ;   < ================================================================================ END noproc=1
; ========================================================================================== END noproc=1

  ss=findfile(star+line_lab+'_resul'+labmac+'.xdr')

; ========================================================================================= BEGIN resul=1
  IF NOT KEYWORD_SET(resul) OR ss(0) EQ '' THEN BEGIN  ;  < =============================== BEGIN resul=1 
; ========================================================================================= BEGIN resul=1

  print,''
  IF ss(0) EQ '' THEN print,star+line_lab+'_resul'+labmac+'.xdr file not found ... lets create it ...'

  ; ------------------------------------------------------------------------------------------------------
  ; --- GOF of the line ----------------------------------------------------------------------------------
  ; ------------------------------------------------------------------------------------------------------

  restore,star+line_lab+'.xdr'

  IF clip EQ 1 THEN BEGIN
     n_ew_down = n_ew_down_cut 
     n_ew_up   = n_ew_up_cut   
  END

  sigma=1./snr

  vs0  = vmin       & Avs  = step_v*vft  & vsf  = fact_vs*vft
  vrt0 = vmin       & Avrt = step_v*vft  & vrtf = fact_vrt*vft

  print,'   (vsmin,vsmax,Avs) = ', vs0,' - ', vsf,' - ', Avs,' km/s',FORMAT='(A23,F5.1,A3,F5.1,A3,F5.1,A5)'
  print,'(vrtmin,vrtmax,Avrt) = ',vrt0,' - ',vrtf,' - ',Avrt,' km/s',FORMAT='(A23,F5.1,A3,F5.1,A3,F5.1,A5)'
  print,''

  ew_in   = ew*(1. - n_ew_down*ew_step)
  ew_fin  = ew*(1. + n_ew_up*ew_step)

  print,'(EW0, EWmin, EWmax, AEW) = ',ew,' - ',ew_in,' - ',ew_fin,' - ',ew_step*ew,' mA',$
         FORMAT='(A27,F5.1,A3,F5.1,A3,F5.1,A3,F5.1,A3)'
  print,''

  vs  = indgen(1+(vsf-vs0)/Avs)*Avs+vs0
  vrt = indgen(1+(vrtf-vrt0)/Avrt)*Avrt+vrt0

  w_inst=w0

  t=w_inst(3)-w_inst(2)

  kk=0

  FOR kew=ew_in,ew_fin,ew_step*ew DO BEGIN
 
    IF kew EQ ew_in THEN print,'USING A DELTA PROFILE AS INTRINSIC PROFILE'

    f_inst=ewgauss(w_inst,kew*1.d-3,200000.d0,xc); Instrumental profile

    IF KEYWORD_SET(em) THEN f_inst=(f_inst-1.)*(-1.)+1.

    PRINT,'EW = ',kew,' (EW_FT0 = ',ew,')',FORMAT='(A5,F5.1,A11,F5.1,A1)'
 
    FOR i=0.d0,n_elements(vs)-1. DO BEGIN
    
      FOR j=0.d0,n_elements(vrt)-1. DO BEGIN

        IF KEYWORD_SET(gauss) THEN BEGIN

          IF j EQ 0 AND i EQ 0 THEN $
	     print,'GOF ... Computing Rotation+Gaussian ... Please, wait for a while ...'
     
          convol,w_inst,f_inst,w_rotrt0,f_rotrt0,vsini=vs(i),res0=0.1*t,vgauss=vrt(j),resol=resol,$
	         nomessage=1,beta=beta

        END ELSE BEGIN

          IF j EQ 0 AND i EQ 0 THEN $
	     print,'GOF ... Computing Rotation+Radial-tangential ... Please, wait for a while ...'

          convol,w_inst,f_inst,w_rotrt0,f_rotrt0,vsini=vs(i),res0=0.1*t,zeta_t=vrt(j),resol=resol,$
	         nomessage=1,beta=beta

        END
     
        f_rotrt = interpol(f_rotrt0,w_rotrt0,wx0)

        chi0a = total(((f_rotrt-fx0)/sigma)^2)/n_elements(wx0)
        chi0 = total((f_rotrt-fx0)^2)
        sigma0 = stddev(f_rotrt-fx0)
        mean0= mean(f_rotrt-fx0)
	
        IF kk EQ 0 THEN sigma_chi = sigma0  ELSE sigma_chi = [sigma_chi,sigma0]
        IF kk EQ 0 THEN chi2      = chi0    ELSE chi2      = [chi2,chi0]
        IF kk EQ 0 THEN chi2_vs   = vs(i)   ELSE chi2_vs   = [chi2_vs,vs(i)]
        IF kk EQ 0 THEN chi2_vrt  = vrt(j)  ELSE chi2_vrt  = [chi2_vrt,vrt(j)]
        IF kk EQ 0 THEN chi2_ew   = kew     ELSE chi2_ew   = [chi2_ew,kew]
      
        kk=kk+1
        
      END

    END

  END

  mm = sort(chi2)

  sigma_chi = sigma_chi(mm)
  chi20     = chi2(mm)
;  chi2      = chi2(mm)*(sigma/min(sigma_chi))^2
  chi2      = chi2(mm)/min(chi2(mm))
  chi2_vs   = chi2_vs(mm)
  chi2_vrt  = chi2_vrt(mm)
  chi2_ew   = chi2_ew(mm)

  save,file=star+line_lab+'_resul'+ labmac+'.xdr',sigma_chi,chi20,chi2,chi2_vs,chi2_vrt,chi2_ew,$
                                                              Avs,vrtf,vsf

; =========================================================================================== END resul=1
  END ;   < ================================================================================= END resul=1
; =========================================================================================== END resul=1

  restore,star+line_lab+'.xdr'
  restore,star+line_lab+'_resul'+labmac+'.xdr'

  w_inst=w0
  t=w_inst(3)-w_inst(2)

  ; ------------------------------------------------------------------------------------------------------
  ; --- Compute mean values and uncertainties from the chi2 distributions --------------------------------
  ; ------------------------------------------------------------------------------------------------------

  chi2=chi2*(min(sigma_chi)*snr)^2    ; TAKE CARE OF THIS LINE !!!!!!!

  final_estimates,chi2_vrt,chi2_vs,chi2_ew,chi2,vft,Avs,vs3,vmft,evmft,vs2,evs2,vg2,evg2,ew2,$
		  chi2_vmac_ft,vsin_ft,vmac_ft,nomessage=1

  env_dist,chi2_vs,chi2,0.,vsf,'vsini',varx_sp_vs,vary_sp_vs,resul_vs
  env_dist,chi2_vrt,chi2,0.,vrtf,'vmac',varx_sp_vm,vary_sp_vm,resul_vm


  ; ------------------------------------------------------------------------------------------------------
  ; --- PLOTS --------------------------------------------------------------------------------------------
  ; ------------------------------------------------------------------------------------------------------

  crea_subplot,[2,2,2],0.10,0.07,xmin=0.1,xmax=0.95,ymin=0.4,pos

  !p.charsize=1.8
  !p.thick=5

  plotin,star+line_lab+labmac+'.ps',mode=1

  ; ======================================================================================================
  ; FIGURE 1: The profile in lambda-space
  ; ======================================================================================================

  ilupal

  ; --- Define the final instrumental profile (w_inst, f_inst) -------------------------------------------

  f_inst=ewgauss(w_inst,ew2*1.d-3,200000.,xc); Instrumental profile

  IF KEYWORD_SET(em) THEN f_inst=(f_inst-1.)*(-1.)+1.

  ; --- Observed profile (maybe renormalized) ------------------------------------------------------------

;  xr=[fix(min(wx0)-0.2*(max(wx0)-min(wx0))),round(max(wx0)+0.2*(max(wx0)-min(wx0))+0.4)]
  xr=[min(wx0)-0.3*(max(wx0)-min(wx0)),max(wx0)+0.3*(max(wx0)-min(wx0))]
  yr=[0.98*min(fx0),1.02*max(fx0)]

  plot,[0],[0],/xsty,/ysty,pos=pos(*,0),title=' ',xtit='!7k!3 [A]',ytit='Norm. flux',yr=yr,xr=xr,xticks=2
  oplot,minmax(w0),1.+(3./snr)*[-1,-1],co=8,line=2
  oplot,w0,f0,co=9,thick=8
  xyerror,wx0,fx0,Ay=fx0*0.+(1./snr)
  simbfill,wx0,fx0,8,0.7,1,fill=1

  ; --- vsini (FT, Theta_RT=0): wp1, fp1 -----------------------------------------------------------------

  convol,w_inst,f_inst,wp1,fp1,vsini=vft,res0=0.1*t,resol=resol,/nomessage,beta=beta
  oplot,wp1,fp1,co=2,line=2

  ; --- vsini (GOF, Theta_RT=0): wp0, fp0 ----------------------------------------------------------------

  convol,w_inst,f_inst,wp0,fp0,vsini=vs3,res0=0.1*t,resol=resol,/nomessage,beta=beta
  oplot,wp0,fp0,co=4

  ; --- vsini (GOF), Theta_RT(GOF): wp2, fp2 -------------------------------------------------------------

  IF KEYWORD_SET(gauss) THEN BEGIN
    convol,w_inst,f_inst,wp2,fp2,vsini=vs2,res0=0.1*t,vgauss=vg2,resol=resol,/nomessage,beta=beta
  END ELSE BEGIN
    convol,w_inst,f_inst,wp2,fp2,vsini=vs2,res0=0.1*t,zeta_t=vg2,resol=resol,/nomessage,beta=beta
  END
  oplot,wp2,fp2,co=3

  ; --- vsini (FT), Theta_RT(GOF, vsini(FT)): wp3, fp3 ---------------------------------------------------

  IF KEYWORD_SET(gauss) THEN BEGIN
    convol,w_inst,f_inst,wp3,fp3,vsini=vft,vgauss=vmft,res0=0.1*t,resol=resol,/nomessage,beta=beta
  END ELSE BEGIN
    convol,w_inst,f_inst,wp3,fp3,vsini=vft,zeta_t=vmft,res0=0.1*t,resol=resol,/nomessage,beta=beta
  END
  oplot,wp3,fp3,co=7,line=2

  ; --- vsini (user), Theta_RT (user): wp4, fp4 ----------------------------------------------------------

  IF KEYWORD_SET(vs_user) THEN BEGIN
    IF vs_user LE 5. THEN vsini=0 ELSE vsini=vs_user
    convol,w_inst,f_inst,wp4,fp4,vsini=vsini,res0=0.1*t,zeta_t=vm_user,resol=resol,/nomessage,beta=beta
    oplot,wp4,fp4,co=6
  END

  ; ======================================================================================================
  ; Summary of info
  ; ======================================================================================================

  IF KEYWORD_SET(gauss) THEN vmac_lab='v!dmac,G!n)' ELSE vmac_lab='v!dmac,RT!n'

  plot,[0],[0],xsty=7,ysty=7,pos=pos(*,2),xtitle=' ',ytitle=' ',yr=[0,1],xr=[0,1]

  xyouts,-0.1,1.10,'Star: '+star,charsize=1.

  xyouts,-0.1,0.95,'Line: '+line+$
	           '   EW='+strtrim(string(ew2,FORMAT='(F5.1)'),2)+' mA'+$
                   '   SNR='+strtrim(string(snr,FORMAT='(I4)'),2),charsize=1.

  xyouts,-0.1,0.70,'FT:  vsini = '+strtrim(string(vft,FORMAT='(F5.1)'),2)+' km/s  ('+vmac_lab+'=0)',$
        charsize=1.0,co=2


  xyouts,-0.1,0.50,'GOF: vsini = '+strtrim(string(vs2,FORMAT='(F5.1)'),2)+' km/s , '+vmac_lab+' = '+$
          strtrim(string(vg2,FORMAT='(F5.1)'),2)+' km/s',charsize=1.0,co=3

  xyouts,-0.1,0.35,'GOF: vsini ('+vmac_lab+'=0) = '+strtrim(string(vs3,FORMAT='(F5.1)'),2)+' km/s',$
          charsize=1.0,co=4

  xyouts,-0.1,0.20,'GOF: '+vmac_lab+' (vsini(FT) = '+strtrim(string(vft,FORMAT='(F5.1)'),2)+') = '+$
          strtrim(string(vmft,FORMAT='(F5.1)'),2)+' km/s',charsize=1.0,co=7

  IF KEYWORD_SET(vs_user) THEN BEGIN
    xyouts,-0.1,0.00,'USER: vsini = '+strtrim(string(vs_user,FORMAT='(F5.1)'),2)+' km/s , '$
          +vmac_lab+' = '+strtrim(string(vm_user,FORMAT='(F5.1)'),2)+' km/s',charsize=1.0,co=6
  END

  xyouts,-0.1, -0.20,'GOF: INTR. PROFILE -> !7d!3 FUNCTION',charsize=1.

  xyouts,-0.1,-0.30,'!7D!3vsini= '+strtrim(string(Avs,FORMAT='(F6.2)'),2)+$
                    ' km/s ,  !7D!3'+vmac_lab+' = '+strtrim(string(Avs,FORMAT='(F6.2)'),2)+' km/s ',$
		    charsize=1.0

  ; ======================================================================================================
  ; FIGURE 2: Fourier transforms of the various profiles
  ; ======================================================================================================

  dd=where(w0 GE min(wx0) AND w0 LE max(wx0))

  wref=w0(dd)

  ; ------------------------------------------------------------------------------------------------------
  ; --- Compute the FT -----------------------------------------------------------------------------------
  ; ------------------------------------------------------------------------------------------------------

  ; --- vsini (FT, Theta_RT=0): wp1, fp1 -----------------------------------------------------------------

  fp1=interpol(fp1,wp1,wref)
  wp1=wref
  
  vfft,wp1,fp1,1.,xc,base,factor,res,nr,max_sigma,ew,sigma_FT0=sigma_FT0

  xfou1=1./(base*factor)
  yfou1=alog10(res(0:nr-1))

  ; --- vsini (GOF), Theta_RT(GOF): wp2, fp2 -------------------------------------------------------------

  fp2=interpol(fp2,wp2,wref)
  wp2=wref

  vfft,wp2,fp2,1.,xc,base,factor,res,nr,max_sigma,ew,sigma_FT0=sigma_FT0

  xfou2=1./(base*factor)
  yfou2=alog10(res(0:nr-1))

  ; --- vsini (FT), Theta_RT(GOF, vsini(FT)): wp3, fp3 ---------------------------------------------------

  fp3=interpol(fp3,wp3,wref)
  wp3=wref

  vfft,wp3,fp3,1.,xc,base,factor,res,nr,max_sigma,ew,sigma_FT0=sigma_FT0

  xfou3=1./(base*factor)
  yfou3=alog10(res(0:nr-1))

  ; --- vsini (user), Theta_RT (user): wp4, fp4 ----------------------------------------------------------

  IF KEYWORD_SET(vs_user) THEN BEGIN

    fp4=interpol(fp4,wp4,wref)
    wp4=wref

    vfft,wp4,fp4,1.,xc,base,factor,res,nr,max_sigma,ew,sigma_FT0=sigma_FT0

    xfou4=1./(base*factor)
    yfou4=alog10(res(0:nr-1))

  END

  vmin=0 & vmax=vsf

  plot,xfou_obs,yfou_obs,xtitle='v sini [km s!u-1!n]',ytitle='ABS (FT AMPLITUDE)',xr=[vmin,vmax],/xstyle,$
       yr=[-4.,0.],/ystyle,title='Fourier transform',pos=pos(*,1)
  oplot,[max_sigma_obs,max_sigma_obs],[-8,-0.7],line=4
  xyouts,max_sigma_obs,-0.6,'Nyq.',charsize=0.7,align=0,orientation=90

  oplot,vft*[1,1],[-20,20],line=2,co=2

  oplot,xfou1,yfou1,co=2,line=2
  oplot,xfou2,yfou2,co=3
  oplot,xfou3,yfou3,co=7,line=3
  IF KEYWORD_SET(vs_user) THEN oplot,xfou4,yfou4,co=6

  ; ======================================================================================================
  ; FIGURE 3: GOODNESS OF FIT (2D)
  ; ======================================================================================================

  !p.title='Goodness-of-fit'

  plot_map,vsf,vrtf,chi2,chi2_vs,chi2_vrt,pos=pos(*,5)

  !p.title=' '

  oplot,[vs2],[vg2],psym=6,co=3
  oplot,[vft],[vmft],psym=6,co=2
  
  IF KEYWORD_SET(vs_user) THEN BEGIN
  
     IF vs_user GT 5. THEN BEGIN

     oplot,[vs_user],[vm_user],psym=6,co=6

     END
     
  END

  ; ======================================================================================================
  ; FIGURE 4: GOF(vsini)
  ; ======================================================================================================

  vmin=0 & vmax=vsf

  !y.range=min(chi2)+[-1.,4.5]
  !x.range=[vmin,vmax]
  !y.style=1

  plot,[0],[0],psym=4,/xsty,/ysty,xtitle='vsini [km s!u-1!n]',ytitle='!7v!3!u2!n',pos=pos(*,3)
  simbfill,chi2_vs,chi2,8,0.8,10,fill=1
  oplot,vs2*[1,1],min(chi2)+[-1,-0.02],co=3,thick=5
  oplot,vft*[1,1],min(chi2)+[-1,-0.02],co=2,line=3,thick=5
  oplot,[vmin,vmax],min(chi2)+[1,1],line=2
  oplot,[vmin,vmax],min(chi2)+[4,4],line=3
  IF KEYWORD_SET(vs_user) AND vs_user GT 5. THEN oplot,vs_user*[1,1],min(chi2)+[-1,-0.02],co=6,line=2,thick=5

  oplot,varx_sp_vs,vary_sp_vs,thick=4,co=0
 
  xyouts,vmin+0.05*(vmax-vmin),min(chi2)-0.4,'!7r!3 (BFM) = '+$
         strtrim(string(min(sigma_chi),FORMAT='(F7.3)'),2),charsize=0.8
  xyouts,vmin+0.05*(vmax-vmin),min(chi2)-0.8,'!7r!3 (SNR) = '+$
         strtrim(string(1./snr,FORMAT='(F7.3)'),2),charsize=0.8

  ; ======================================================================================================
  ; FIGURE 5: GOF(Theta_RT)
  ; ======================================================================================================

  tt=where(chi2 LE min(chi2)+16)
 
  vmin=0 & vmax=max(chi2_vrt(tt))+min(chi2_vrt(tt))
  
  !x.range=min(chi2)+[4.5,-1]
  !y.range=[vmin,vmax]

  IF KEYWORD_SET(gauss) THEN ytit='v!dmac!n (G) [km s!u-1!n]' ELSE ytit='v!dmac!n (R-T) [km s!u-1!n]'

  plot,[0],[0],/nodata,/xsty,/ysty,ytitle=ytit,xtitle='!7v!3!u2!n',pos=pos(*,4)
  simbfill,chi2,chi2_vrt,8,0.8,10,fill=1
  oplot,min(chi2)+[-1,-0.02],vg2*[1,1],co=3,thick=5
  oplot,min(chi2)+[-1,-0.02],vmft*[1,1],co=7,line=3,thick=5
  oplot,min(chi2)+[1,1],minmax(chi2_vrt),line=2
  oplot,min(chi2)+[4,4],minmax(chi2_vrt),line=3
  IF KEYWORD_SET(vs_user) AND vs_user GT 5. THEN oplot,min(chi2)+[-1,-0.02],vm_user*[1,1],co=6,line=2,thick=5

  oplot,vary_sp_vm,varx_sp_vm,thick=4,co=0

  plotout

  IF KEYWORD_SET(ps) THEN spawn, ps_viewer+' '+star+line_lab+labmac+'.ps&'

  ; ------------------------------------------------------------------------------------------------------
  ; --- The final results are printed in the terminal ----------------------------------------------------
  ; ------------------------------------------------------------------------------------------------------

  print,''
  print,'====================================== FINAL RESULTS ============================================'
  print,''

  print,'STAR: '+star
  print,''
  print,'INTRINSIC PROFILE: DELTA FUNCTION'+'  LINE: '+line

  print,''
  print,'EW (mA) SNR  vsini(FWHM)  vsini(FT)   Theta_RT(FT+GOF)       vsini(GOF)         Theta_RT(GOF)'
  print,'------- ---  -----------  ---------  ------------------   ------------------   ------------------'
  print,string(ew2,FORMAT='(I4)')+'   '+$
        string(snr,FORMAT='(I4)')+'     '+$
        string(vs3,FORMAT='(F5.1)')+'      '+$
        string(vft,FORMAT='(F5.1)')+'    '+$
        string(vmft,FORMAT='(F5.1)')+' - '+string(evmft(0),FORMAT='(F4.1)')+$
                                     ' + '+string(evmft(1),FORMAT='(F4.1)')+'  '+$
        string(vs2,FORMAT='(F5.1)')+' - '+string(evs2(0),FORMAT='(F4.1)')+$
                                    ' + '+string(evs2(1),FORMAT='(F4.1)')+'  '+$
        string(vg2,FORMAT='(F5.1)')+' - '+string(evg2(0),FORMAT='(F4.1)')+$
                                    ' + '+string(evg2(1),FORMAT='(F4.1)')
  print,''
  print,'====================================== FINAL RESULTS ============================================'
  print,''

  d={star: star, line:line, ew:ew2, snr:snr, vs0: vs3, vft:vft, vmft:vmft, evmft:evmft, vsgof:vs2, $
     evsgof: evs2, vmgof:vg2, evmgof:evg2}

  ; -------------------------------------------------------------------------------------------------------

  !p.thick=1

end

; =========================================================================================================
; =========================================================================================================
;                                  SOME AUXILIARY PROGRAMS 
; =========================================================================================================
; =========================================================================================================

  pro smooth_edge,wx,fx,bb,perc=perc
  
    fx  = fx-bb
    nx  = n_elements(wx)

    nxi = nx       ; Number of sampling points
   
    win_cos,nx,correc,perc=perc

    fx   =   fx*correc
    fx   =   fx+bb

  end

; =========================================================================================================
; =========================================================================================================

  pro extend_cont,wx,fx,ampli,bb
  
    nx  = n_elements(wx)

    plus = ampli*nx

    IF plus ne 0 THEN BEGIN

      paso      =  wx(1)-wx(0)

      wextens1  =  (findgen(plus)+1.)*paso+wx(nx-1)
      wextens2  =  -(-findgen(plus)+plus)*paso+wx(0)

      wx0       =  wx
      wx        =  [wextens2,wx,wextens1]

      fextens1  =  fltarr(plus)+bb
      fextens2  =  fltarr(plus)+bb

      fx0       =  fx
      fx        =  [fextens1,fx,fextens2]

    END

  end

; =========================================================================================================
; =========================================================================================================

  pro win_cos,nx,wx,perc=perc,help=help

    wx=replicate(1.,nx)
   
    IF KEYWORD_SET(perc) THEN BEGIN
   
      nxw=nint(perc*nx/100.)
   
      wi=0.5*(1.-cos(!pi*findgen(nxw)/nxw))
      wx(0:nxw-1)=wi
      wx(nx-nxw:nx-1)=reverse(wi(0:nxw-1))

    END
      
  end

; =========================================================================================================
; =========================================================================================================

  pro vfft,wx,fx,ycont,xc,base,factor,res,nr,max_sigma,ew,nyq=nyq,sigma_FT0=sigma_FT0

    IF NOT KEYWORD_SET(sigma_FT0) THEN sigma_FT0=0.660d0

    ycont=1

    perc   =         10   ; % of the line sampling that is smoothed
    ampli  =         40   ; Extension factor for the line continuum
    cc     =     2.99e5   ; Light speed

    wf=wx
    ff=fx

    smooth_edge,wf,ff,ycont,perc=perc
  
    extend_cont,wf,ff,ampli,ycont


    nx   =  n_elements(wf)
    t    = wf(3)-wf(2)
    ds   = 1.d0/(t*nx)
    max_sigma=0.5d0/t

    res = abs(fft(ff-ycont))
    ew  = (res(0)*nx*t)
    nr= n_elements(res)-1.

    res=res/res(0)

    factor=xc/cc/sigma_FT0

    base=[0.,(indgen(nr-1)+1)*ds]

    xrmax=0.8/t
    xrmin=base(min(where(res(0:nr-1) LT 0.9)))

    xrmax=1./(xrmax*factor)
    xrmin=1./(xrmin*factor)
    max_sigma=1./(max_sigma*factor)

    IF KEYWORD_SET(nyq) THEN print,'Minimum vsini (Nyquist): ', max_sigma, ' km/s',FORMAT='(A25,F5.1,A5)'

  end

; =========================================================================================================
; =========================================================================================================

  pro plot_lin,lam_list

    FOR i=0,n_elements(lam_list)-1 DO BEGIN

      oplot,lam_list(i)*[1.,1.],[0,10],line=2

    END

  end

; =========================================================================================================
; =========================================================================================================

  pro final_estimates,chi2_vrt,chi2_vs,chi2_ew,chi2,vft,Avs,vs3,vmft,evmft,vs2,evs2,vg2,evg2,ew2,$
	   	      chi2_vmac_ft,vsin_ft,vmac_ft,nomessage=nomessage

    IF NOT KEYWORD_SET(nomessage) THEN BEGIN

      print,''
      print,' vsini(FT) = ',vft,FORMAT='(A28,F5.1)'

    END
  
  ; --- vsini (GOF , Theta_RT=0) --- ** vs3 **

    mm=where(chi2_vrt EQ min(chi2_vrt))
    vfwhm=chi2_vs(mm)

    chi2_vfwhm=chi2(mm)  ; chi2 distribution for Theta_RT=0

    env_dist,vfwhm,chi2_vfwhm,min(vfwhm),max(vfwhm),'vsini(FWHM)',varx_sp_vs,vary_sp_vs,resul

    vs3=resul.bfm_1x

    IF NOT KEYWORD_SET(nomessage) THEN BEGIN
  
      print,''
      print,'vsini(GOF, Theta_RT=0) = ',vs3,FORMAT='(A28,F5.1)'

    END

  ; --- Theta_RT (GOF , vsini=vft)  ** vmft **

    mm=where(chi2_vs GE vft-0.5*Avs AND chi2_vs LE vft+0.5*Avs)
    vsin_ft=chi2_vs(mm)
    vmac_ft=chi2_vrt(mm)

    chi2_vmac_ft=chi2(mm)  ; chi2 distribution for vsini=vsini(FT)

    env_dist,vmac_ft,chi2_vmac_ft,min(vmac_ft),max(vmac_ft),'vmac(vsini,FT)',varx_sp_vs,vary_sp_vs,resul

    vmft=resul.bfm_1x
    evmft=[resul.errmin_1x,resul.errmax_1x]
    
    IF NOT KEYWORD_SET(nomessage) THEN BEGIN

      print,'Theta_RT (GOF, vsini(FT)) = ',vmft,' - ',evmft(0),' + ',evmft(1),FORMAT='(A28,F5.1,2(A3,F5.1))'

    END

  ; --- vsini, Theta_RT (GOF) ** vs2, vg2 **

    env_dist,chi2_vs,chi2,min(chi2_vs),max(chi2_vs),'vsini (GOF)',varx_sp_vs,vary_sp_vs,resul

    vs2=resul.bfm_1x
    evs2=[resul.errmin_1x,resul.errmax_1x]
    
    env_dist,chi2_vrt,chi2,min(chi2_vrt),max(chi2_vrt),'vmac (GOF)',varx_sp_vm,vary_sp_vm,resul

    vg2=resul.bfm_1x
    evg2=[resul.errmin_1x,resul.errmax_1x]
    
    env_dist,chi2_ew,chi2,min(chi2_ew),max(chi2_ew),'EW (GOF)',varx_sp_vs,vary_sp_vs,resul

    ew2=resul.bfm_1x

    IF NOT KEYWORD_SET(nomessage) THEN BEGIN

      print,' vsini (GOF) = ',vs2,' - ',evs2(0),' + ',evs2(1),FORMAT='(A28,F5.1,2(A3,F5.1))'   
      print,' Theta_RT (GOF) = ',vg2,' - ',evg2(0),' + ',evg2(1),FORMAT='(A28,F5.1,2(A3,F5.1))'
      print,' EW (GOF) = ',ew2,FORMAT='(A28,F5.1)'      
      print,''
  
    END

end

; =========================================================================================================
; =========================================================================================================

pro plot_map,vsf,vrtf,chi2,chi2_vs,chi2_vrt,pos=pos,gauss=gauss

  chi0=chi2-min(chi2) 

  ilupal

  tt=where(chi2 LE min(chi2)+9)
  
  vmin=0 & vmax=vsf                & !x.range=[vmin-0.05*(vmax-vmin),1.02*vmax]
  vmin=0 & vmax=max(chi2_vrt(tt))+min(chi2_vrt(tt))  & !y.range=[vmin-0.05*(vmax-vmin),1.02*vmax]

  IF KEYWORD_SET(gauss) THEN ytit='v!dmac!n (G) [km s!u-1!n]' ELSE ytit='v!dmac!n (R-T) [km s!u-1!n]'
  
  plot,[0],[0],/nodata,/xsty,/ysty,xtitle='v sini [km s!u-1!n]',ytitle=ytit,pos=pos

  contour,chi0,chi2_vs,chi2_vrt,/irregular,/overplot ,levels=min(chi0)+[0,1,4,9],cell_fill=1,$
  c_colors=[9,10,11,1]

  contour,chi0,chi2_vs,chi2_vrt,/irregular,/overplot ,levels=min(chi0)+[0,1,4,9],$
  c_labels=[0,1,1,1],c_annotation=['','1!7r!3','2!7r!3','3!7r!3'],c_charsize=1.2,c_charthick=3

  oplot,min(chi2_vs)*[1,1],[0,vrtf],line=2
  oplot,[0,vsf],min(chi2_vrt)*[1,1],line=2
  
end

; =========================================================================================================
; =========================================================================================================

  pro env_dist,var,chi2,vmin,vmax,title,varx_sp,vary_sp,resul,inv=inv

  ; -------------------------------------------------------------------------------------------------------
  ; Curvature for the spline fitting (high values means linear fit)
  ; -------------------------------------------------------------------------------------------------------

  curv=5.   ; 5.

  ; -------------------------------------------------------------------------------------------------------
  ;  Struture in which the results for a given variable will be stored
  ; -------------------------------------------------------------------------------------------------------

  resul={title:'' ,bfm:0. ,min_1x :0. , max_1x :0. , bfm_1x :0. ,errmin_1x :0. ,errmax_1x :0. ,err_1x :0. ,$
                           min_2x :0. , max_2x :0. , bfm_2x :0. ,errmin_2x :0. ,errmax_2x :0. ,err_2x :0. ,$
	  	 	   var_min:0. , var_max:0. }

  resul.title = title

  ; -------------------------------------------------------------------------------------------------------
  ;  Some info for the plots
  ; -------------------------------------------------------------------------------------------------------

  cs1 = 0.7 & cs2 = 1 & cs3 = 1 & csp=1.7

  var=double(var)

  var_min=vmin       & var_max=vmax        ; Minimum and maximum values for the plots (x-axis)
  ymin=min(chi2)-1.  & ymax=min(chi2)+5.   ; Minimum and maximum values for the plots (y-axis)

  resul.var_max = var_max
  resul.var_min = var_min
 
  xsize=var_max-var_min

  nstep=1000. & sh=(var_max-var_min)/nstep ; Quantities used for the splines fitting 

  x1=var & y1=chi2 
  mm=sort(x1) & x1=x1(mm) & y1=y1(mm)

  ; -------------------------------------------------------------------------------------------------------
  ; Fit to the envelope of the distribution + uncertainty estimation
  ; -------------------------------------------------------------------------------------------------------

  var_x = var[UNIQ(var, SORT(var))]
  var_y = var_x

  FOR in=0,n_elements(var_x)-1 DO BEGIN

    mm=where(var EQ var_x(in))
    var_y(in)=min(chi2(mm))

  END

  ff0=where(chi2 LE min(chi2))

  resul.bfm=var(ff0)

;  print,resul.title+' [BFM] : ',resul.bfm,FORMAT='(A,A,F5.2,2X)'

  IF n_elements(var_x) EQ 1 THEN resul.bfm_1x    = var_x(0)
  IF n_elements(var_x) EQ 1 THEN resul.bfm_2x    = var_x(0)

  IF n_elements(var_x) GE 3  THEN BEGIN

    varx_sp=var_min+findgen(nstep)*sh
    vary_sp=spline(var_x,var_y,varx_sp,curv)

    ; -----------------------------------------------------------------------------------------------------
    ; -------- RESULTS FOR 1sigma -----------------
    ; -----------------------------------------------------------------------------------------------------

    ff1=where(vary_sp LE min(y1)+1.)

    IF ff1(0) GE 0 THEN BEGIN
  
    resul.min_1x    = min(varx_sp(where(vary_sp LE min(vary_sp)+1.)))   ;  min(varx_sp(ff1))
    resul.max_1x    = max(varx_sp(where(vary_sp LE min(vary_sp)+1.)))   ;  max(varx_sp(ff1))
    resul.bfm_1x    = varx_sp(where(vary_sp LE min(vary_sp))) ; 0.5*total(minmax(varx_sp(ff1)))
    resul.errmin_1x = abs(resul.min_1x-resul.bfm_1x)
    resul.errmax_1x = abs(resul.max_1x-resul.bfm_1x)
    resul.err_1x    = abs(resul.min_1x-resul.bfm_1x)

    END
  
  ; -------------------------------------------------------------------------------------------------------
  ; -------- RESULTS FOR 2sigma -----------------
  ; -------------------------------------------------------------------------------------------------------

    ff2=where(vary_sp LE min(y1)+4.)

    IF ff2(0) GE 0 THEN BEGIN
  
    resul.min_2x    = min(varx_sp(ff2))
    resul.max_2x    = max(varx_sp(ff2))
    resul.bfm_2x    = varx_sp(where(vary_sp LE min(vary_sp))) ; 0.5*total(minmax(varx_sp(ff2)))
    resul.errmin_2x = abs(resul.min_2x-resul.bfm)
    resul.errmax_2x = abs(resul.max_2x-resul.bfm)
    resul.err_2x    = abs(resul.min_2x-resul.bfm_2x)

    END
  
  END

;  print,'min, max, mean [1sigma]:', resul.min_1x,resul.max_1x,resul.bfm_1x,FORMAT='(A,3(F7.2,2X))'

end
