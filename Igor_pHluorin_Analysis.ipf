#pragma rtGlobals=1		

Function LoadStack_TILL(stackname, pathname, filename,frameN, timewave,timestep) 


 	String stackname, pathname, filename
 	Variable frameN
	String timewave
	Variable timestep 
 
 
 
 	Variable sti_file, i, pnts, h, m, s, ss, day, month, year, sec_offset, ss_offset
 	String str, buf1, buf2
 	Variable /G background_buffer=0 												
 
 	
	NewPath /O/Q sti_path, pathname 												
																				


	String extension 
	Variable chars=4, start_frame=0, end_frame=frameN-1, curr_image                 

	Make /O/N=(frameN) /T $stackname											

	Wave /T names=$stackname



	for (curr_image=start_frame;curr_image<=end_frame;curr_image+=1)						

			extension="_"+num2str(curr_image)									

			names[curr_image-start_frame]=stackname+extension					

			ImageLoad /Q/T=tiff /O/P=sti_path /N=$(names[curr_image-start_frame])  filename+extension+".tif"  

	endfor
	
	


	Variable k=DimSize($(names[0]),0)												
	Variable l,j												

	printf "Scaling images..."


	if (mod(k,2)!=0)																
	
  			for(i=0;i<frameN;i+=1)
  			
   				DeletePoints /M=0 k-1, 1, $(names[i])								
   				
  			endfor
	endif 


	k=DimSize($(names[0]),1)	
																				

	if (mod(k,2)!=0)
	
 			for(i=0;i<frameN;i+=1)
 			
  			 DeletePoints /M=1 k-1, 1, $(names[i])									
  			 
 			 endfor
	endif 



FindScale(stackname)															

	NVAR min_scale=min_imscale  												

	NVAR max_scale=max_imscale  

	printf "...done\r"																



Display as stackname+" (1/"+num2str(frameN)+")";AppendImage $(names[0])			

ModifyImage $(names[0]) ctab= {(min_scale),(max_scale),Grays,0}						

SetAxis/A/R left																	

DoWindow/C $(stackname+"_win")													


KillPath sti_path 			

Make /O/N=(frameN) $timewave													
													
Wave xwave=$timewave

xwave=p*timestep																
 
End

//
//
//
//
//


Function FindScale(stackname)



String stackname


NVAR min_scale=min_imscale  



     	if (NVAR_exists(min_scale)==0)
     	
       	Variable /G min_imscale
       	
     		NVAR min_scale=min_imscale
     		
     endif 

min_scale=2^32



NVAR max_scale=max_imscale  

     	if (NVAR_exists(max_scale)==0)
     	

     		
      	 	Variable /G max_imscale
       
     		NVAR max_scale=max_imscale
     	   	
      endif 

max_scale=0



Wave /T names=$stackname

Variable frameN=numpnts(names),i,j,l



for(i=0;i<frameN;i+=1)

	Wave imm=$(names[i])
	
	WaveStats /Q imm
 
 	if(min_scale>V_min)
   		min_scale=V_min
 	endif
 
 	if(max_scale<V_max)
   		max_scale=V_max
 	endif
 	
endfor



printf "\r Minimal stack intensity: %d", min_scale

printf "\r Maximal stack intensity: %d", max_scale


printf "\r\r"

End

//
//
//
//
//


Function subtractvalfromwvs(wvpath,wvnmebase,wvnum,bgval)					


string wvpath,wvnmebase

variable wvnum, bgval



variable i



for (i=0;i<wvnum;i+=1)
	
	wave wv0=$wvpath+wvnmebase+num2str(i)
	
	redimension/D wv0													
	
	wv0-=bgval
	
endfor

end


//
//
//
//
//


Function kinetics_maskonly(stackname,baslst, baslen, maxist, maxien)		



String stackname

Variable baslst, baslen, maxist, maxien  



Wave/T stk=$stackname																		

Variable frameN=numpnts(stk),i											

AverageStack(stackname, baslst+1, baslen+1)								

Duplicate /O $("aver_"+stackname) basl_imag								

AverageStack(stackname, maxist+1, maxien+1)								

Duplicate /O $("aver_"+stackname) maxi_imag								

Duplicate /O maxi_imag testtemplate 										

testtemplate=maxi_imag-basl_imag											

ROI_stack (5, 3, "testtemplate", stackname, 15,100,"timewave") 			

End


//
//
//
//
//


Function AverageStack(stackname, st, en) 									



String stackname

Variable st, en 															



Wave/T stk=$stackname													

Wave im=$(stk[0])														

Variable N=(en-st)+1, i,y=DimSize(im,0), x=DimSize(im,1)					

if (N==0)	
																
 		DoAlert 0, "No images to average!"									
 
 		return -1
endif

Make /D/O/N=(y,x) $("aver_"+stackname)									

Wave aver_im=$("aver_"+stackname)
 
CopyScales im, aver_im													

aver_im=0																

for (i=st;i<=en;i+=1)														

		 Wave im=$(stk[i-1])												
		 
 		aver_im+=im/N													
endfor

End

//
//
//
//
//

Function ROI_stack(waveletdepth, waveletstart, template, stackname, arealowthresh,areahighthresh,timing) 



Variable waveletdepth, waveletstart

String template, stackname

Variable arealowthresh, areahighthresh

String timing




Wave img=$template																	

Wave /T stk=$stackname													
String s=template+"_mask"												

staWavelet3(waveletdepth, waveletstart,img, s)								


Wave msk=$s															

Variable frameN=numpnts(stk), i, j, k, l, dimx, dimy, sx,sy1, sy2, xtemp, ytemp	

Make /O/N=(1) part_num

Make /O/T/N=(1) part_coordX, part_coordY, part_area, part_circ, part_index, part_intens1,part_intens2

print "Done Wavelet"

i=0

part_coordX[i]=template+"_partX"											

part_coordY[i]=template+"_partY"

part_area[i]=template+"_partarea"

part_circ[i]=template+"_partcirc"

part_index[i]=template+"_partindex"

ImageAnalyzeParticles /A=(arealowthresh)/Q stats msk						

part_num[i]=V_NumParticles												
Duplicate/O W_circularity $(part_circ[i])									

Duplicate/O W_ImageObjArea $(part_area[i])

Duplicate/O/D W_SpotX $(part_coordX[i])

Duplicate/O/D W_SpotY $(part_coordY[i])

Make /O/W/N=(V_NumParticles) $(part_index[i])								

Make /O/D/T/N=(V_NumParticles) $("intens_index")							

Wave index=$(part_index[i])												
							
Wave /T intens=$("intens_index")											

index=-1
  
  Wave p_X=$(part_coordX[i]), p_Y=$(part_coordY[i]), p_area=$(part_area[i]), p_circ=$(part_circ[i]) 
   
 for(j=0; j<numpnts(intens);j+=1)											
 
 intens[j]="ROI_"+num2str(j)												
 
 Make/N=(frameN)/O/D $(intens[j])
 
 Wave timetrace=$(intens[j])												
 
 if (p_area[j]>areahighthresh)												
 
   KillWaves /Z $(intens[j])
   
   DeletePoints j,1, intens, p_area,p_circ,p_X,p_Y,index
   
   j-=1
   continue
   
   
 endif
 xtemp=p_X[j]; ytemp=p_Y[j] 
 ImageAnalyzeParticles /L=(xtemp,ytemp)/Q mark msk
 
 Duplicate /O/D M_ParticleMarker temp_mask
 temp_mask=1-temp_mask/64
 
   Redimension/B/U temp_mask
ImageMorphology/O/E=5 Dilation temp_mask
temp_mask/=2

 
 for(i=0;i<frameN;i+=1)
   Wave img=$(stk[i])
  dimx=DimSize(temp_mask,0); dimy=DimSize(temp_mask,1)
  sx=0;sy1=0;sy2=0
   for(k=max(0,xtemp-30);k<min(dimx,xtemp+30);k+=1)
   for(l=max(0,ytemp-30);l<min(dimy,ytemp+30);l+=1)
     sy1+=img[k][l]*temp_mask[k][l]
   endfor
   endfor
    timetrace[i]=sy1/p_area[j]
    
 endfor

endfor


KillWaves /Z temp_mask


Display;AppendImage msk

Edit p_area, p_circ, p_X, p_Y

InsertPoints 0,1, intens
intens[0]=timing



End

//
//
//
//
//

Function staWavelet3(depth, start,image, outmask)



variable depth, start

wave image															

String outmask															



	
Duplicate /D/O image imgstack											      
	
silent 1																	
	
variable nchar, index=0

variable cnt=0, maxdepth, nth=0, ind1, ind2, void=0
	
variable lenx = dimsize(imgstack,0)												
variable leny = dimsize(imgstack,1)											
								
maxdepth =trunc( ( log( (min(lenx,leny)-3)/2+1) / log(2) )+1)					
	
	
	if (depth>maxdepth)												
												
		depth=maxdepth
		
		print depth
		
	endif


	
variable ngaps, newdim, Ithresh

variable ii=1, jj=0, mlength,nn,ss
	
make/D/O/N=(lenx, leny) imgAii=0, imgWt=0, $(outmask), imgfiltered=0			

make/D/O/N=(lenx, leny, depth+1) imgA=0									

make/D/O/N=(lenx, leny, depth) imgW=0, imgWthard=0						

make/D/O/N=(lenx, leny, depth-start) imgP=1								

make/D/O/N=3 kernel={1/8,3/4,1/8}											

Wave imgmask=$(outmask)												


	
ngaps=0

newdim=0

Ithresh=0

ii=1

jj=0

mlength=0

nn=0

ss=0	



	imgA[][][0] = imgstack[p][q][nth]										
	
	imgfiltered[][] = imgstack[p][q][nth]										
	
	do
		ngaps=(2^(ii-1)-1)																			
		
		newdim=3+2*ngaps												
		
		duplicate/O kernel, kernelj											
		
		insertpoints 2, ngaps, kernelj										
		
		insertpoints 1, ngaps, kernelj										
	
		make/D/O/N=(newdim, newdim) filter=1								
		
		filter*=kernelj[p]*kernelj[q]											
		
		imgAii = imgA[p][q][ii-1]											
		
		matrixconvolve filter imgAii											
		
		imgA[][][ii] = imgAii[p][q]											
		
		imgW[][][ii-1]=imgA[p][q][ii-1]-imgA[p][q][ii]							
		
		imgWt = imgW[p][q][ii-1]											
		
		void = Wthreshold(imgWt)										
		
		imgWthard[][][ii-1]=imgWt[p][q]										
		
		if (ii>start)														
		
			imgP[][][ii-start-1]=imgP[p][q][ii-start-2]*imgWt[p][q]			
			
		endif
		
		ii+=1
		
	while(ii<(depth+1))

	imgmask=imgP[p][q][depth-start-1]										
	
	duplicate/O/D imgmask Maskstats										
	
	Mlength=dimsize(imgmask,0)*dimsize(imgmask,1)						
	
	redimension/N=(Mlength) Maskstats									
	
	wavestats/Q Maskstats												
	
	Maskstats=abs(V_avg-Maskstats[p])									
							
	sort Wstats, Wstats													
	
	Ithresh = 3*Maskstats[mod(Mlength, 2)]/0.67								
	
		void = imgbinary(imgmask, Ithresh)									
		
		void = removepxls(imgmask, 2)										
		
		void = removepxls(imgmask, 2)										
		
		void = removepxls(imgmask, 2)									
	
	imgmask-=1															
	
	imgmask*=-255														
	
	Redimension/B/U imgmask											

End

//
//
//
//
//

function/D Wthreshold(Wimg)


wave/D Wimg



duplicate/O/D Wimg Wstats											

variable length=dimsize(Wimg,0)*dimsize(Wimg,1)						

	redimension/N=(length) Wstats										
	
	wavestats/Q Wstats											
	
	Wstats=abs(V_avg-Wstats[p])									
	
	sort Wstats, Wstats												
	
	variable thard = 3*Wstats[mod(length, 2)]/0.67						
	
	variable/D xx=0, yy
	
	do										
		yy=0
		
		do
			if (Wimg[xx][yy]<thard)									
			
				Wimg[xx][yy]=0
				
			endif
			
			yy+=1													
			
		while(yy< dimsize(Wimg,1))
		
		xx+=1														
		
	while(xx< dimsize(Wimg,0))
	
end

//
//
//
//
//

function/D imgbinary(Pimg, thresh)



wave/D Pimg

variable thresh




variable/D xx=0, yy

	do
		yy=0
		
		do
			if (Pimg[xx][yy]<thresh)
			
				Pimg[xx][yy]=0									
				
			else
			
				Pimg[xx][yy]=1				
			endif
				
			yy+=1	
								
		while(yy< dimsize(Pimg,1))
		
		xx+=1	
		
	while(xx< dimsize(Pimg,0))
end

//
//
//
//
//

function/D removepxls(Pimg, numneighb)



wave/D Pimg

variable numneighb													



make/O/N=(3,3) Anb=0

variable/D xx=0, yy, Anbsum
	do
		yy=0
		do
			Anb[][]=Pimg[xx-1+p][yy-1+q]
			Anbsum=sum(Anb,0,8)
			 
			if (Anbsum<=(numneighb+1))
				Pimg[xx][yy]=0
			endif
			yy+=1
		while(yy< dimsize(Pimg,1))
		xx+=1	
	while(xx< dimsize(Pimg,0))
end
//
//
//
//
//

Function AverageWaves(index_name, aver_name)



String index_name, aver_name
 


 Wave /T index=$index_name
 
 Wave xwav=$(index[0])
 
 Variable i,j, M=numpnts(index)-1, N=numpnts(xwav)
 
 Make /O/N=(N) $aver_name
 
 Wave average=$aver_name
 
 Duplicate /O average, $(aver_name+"_stdev")
 
 Wave stdev=$(aver_name+"_stdev")
 
 Make /O/N=(M) buffer
 
 for (i=0;i<N;i+=1)
 
   	for (j=1;j<=M;j+=1)
   	
    		Wave w=$(index[j])
    		
    		buffer[j-1]=w[i]
   	endfor
   	
  	Wavestats /Q buffer
  
  	average[i]=V_avg
  
  	stdev[i]=V_sdev
  	
 endfor
 
 DoWindow /K $(aver_name+"_graph")
 
 Display average vs xwav as "Graph of "+aver_name
 
 ErrorBars $aver_name,Y wave=(stdev,stdev)
 
 DoWindow /C $(aver_name+"_graph")
 
 End
 //
 //
 //
 //
 //
 
Function NormWaves(index_name, averwave)



 String index_name, averwave


 
 Wave /T index=$index_name
 
 Wave xwav=$(index[0])
 
 Wave average=$averwave
 
 Variable x0, x1
 


 DoWindow /K $"tmpwindow"
 
 Display average vs xwav as "Select the baseline with cursors"
 
 DoWindow /C $"tmpwindow"
 
 ShowInfo
 
 DummyPanel()
 
 x0=pcsr(A)
 
 x1=pcsr(B)
 
 Variable j, M=numpnts(index)-1, N=numpnts(xwav)
 
 Duplicate /O xwav $(index[0]+"_norm")
 
 Make /O/T/N=(M+1) $(index_name+"_norm")
 
 Wave /T index_norm=$(index_name+"_norm")
 
 index_norm[0]=index[0]+"_norm"
 
   for (j=1;j<=M;j+=1)
   
    		Duplicate /O $(index[j]) $(index[j]+"_norm")
    		
    		Wave w=$(index[j]+"_norm")
    	
    	 	Wavestats /Q/R=[x0,x1] w
    	 
    	 	w/=V_avg
    	 
    		index_norm[j]=index[j]+"_norm"
    	
   endfor
 
 DoWindow /K $"tmpwindow"
  
 Edit index_norm
 
 End
 //
 //
 //
 //
 //
 
 Function NormWaves_NH4(index_name, averwave)



 String index_name, averwave
 
 
 
 Wave /T index=$index_name
 
 Wave xwav=$(index[0])
 
 Wave average=$averwave
 
 Variable x0, x1, x2, x3
 
 


 DoWindow /K $"tmpwindow"
 
 Display average vs xwav as "Select the baseline with cursors"
 
 DoWindow /C $"tmpwindow"
 
 ShowInfo
 
 DummyPanel()
 
 x0=pcsr(A)
 
 x1=pcsr(B)
 
 DoWindow /K $"tmpwindow"
 
 Display average vs xwav as "Select the NH4 baseline for subtraction with cursors"
 
 DoWindow /C $"tmpwindow"
 
 ShowInfo
 
 DummyPanel()
 
 x2=pcsr(A)
 
 x3=pcsr(B)
 
 Variable j, M=numpnts(index)-1, N=numpnts(xwav)
 
 Duplicate /O xwav $(index[0]+"_norm")
 
 Make /O/T/N=(M+1) $(index_name+"_norm")
 
 Wave /T index_norm=$(index_name+"_norm")
 
 index_norm[0]=index[0]+"_norm"
 
   for (j=1;j<=M;j+=1)
   
    	Duplicate /O $(index[j]) $(index[j]+"_norm")
    	
    	Wave w=$(index[j]+"_norm")
    	
     	Wavestats /Q/R=[x2,x3] w
     	
     	w-=V_avg
     	
     	Wavestats /Q/R=[x0,x1] w
     	
    	 w/=V_avg
    	 
  	 index_norm[j]=index[j]+"_norm"
  	 
   endfor
 
DoWindow /K $"tmpwindow"
  
Edit index_norm
 
End
 //
 //
 //
 //
 //
Function SelectWaves2(list_name,selected1, selected2)



      Wave/T list_name
      
      String selected1, selected2
      

      
    Variable i,Nsel1, Nsel2
      
      Wave xwav=$(list_name[0])
      
      Make/O/N=1/T $selected1, $selected2
      
      Wave /T sel1=$selected1, sel2=$selected2
      
      sel1[0]=list_name[0]
      
      sel2[0]=list_name[0]
          
       Wave /T reg_done=$list_name(1)
          
       for (i=1;i<numpnts(list_name);i+=1)
      
      			Wave inw=$list_name[i]
      
        		Nsel1=numpnts(sel1)
        
        		Nsel2=numpnts(sel2)
        
    			DoWindow/K tmpwindow
    
        		Display inw vs xwav
        		
        		ModifyGraph nticks(bottom)=30
        		
        		DoWindow/C tmpwindow
        		
        		ModifyGraph nticks(bottom)=30
       
        		DoUpdate    
        		
    			DoAlert 2, list_name[i]+"looks synaptic (->yes), non-synaptic (->no), no release (->cancel)"
    
    			If (V_Flag==1)
    
        			InsertPoints Nsel1,1, sel1
        
       			 sel1[Nsel1]=list_name[i]
        
    			elseif (V_Flag==2)
    
       			 InsertPoints Nsel2,1, sel2
        
        			sel2[Nsel2]=list_name[i]
        
    			elseif (V_flag==3)
    
       		continue
       		
      			endif
      
     	endfor

DoWindow/K tmpwindow

ModifyGraph nticks(bottom)=30

printf "%d time traces are synaptic out of %d\r", numpnts(sel1)-1, numpnts(list_name)-1

printf "%d time traces are detonators out of %d\r", numpnts(sel2)-1, numpnts(list_name)-1
        
End
//
//
//
//
//
Function DummyPanel()

      NVAR drawing_mode=$"drawings_colors"
      
      
	DoWindow/F Buttons
	
	if( V_Flag==1 )
	
		return 0
	endif
     
      	NewPanel /K=1 /W=(502.8,417.2,693,420.8) as "Kill me when done!!!"
      	
	DoWindow/C Buttons
	
	AutoPositionWindow/E/M=1

	ModifyPanel fixedsize=1
	
	PauseForUser Buttons, $"tmpwindow"

end

