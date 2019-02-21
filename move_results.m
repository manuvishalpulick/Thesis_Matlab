function move_results(mk,Tmp)



    if Tmp == 0.000
           movefile('Evolution_subplot*',mk)                 % move all the video files to that directory
          %movefile('filmEvolution*',mk)                 % move all the video files to that directory
          %movefile('Fourier_evolution*',mk)                 % move all the video files to that directory
           movefile('harmonic_evolution*',mk)                 % move all the image files to that directory
           movefile('Final_Sk_snapshot*',mk)                 % move all the image files to that directory
           movefile('Maximum_energy_frequency*',mk)                 % move all the image files to that directory
           movefile('Rupture_instant*',mk)                 % move all the image files to that directory 
%           movefile(strcat('*_',num2str(L_flat),'_*.mat'),mk)                 % move the .mat file to that directory
%           movefile('wave_num_fig*',mk)                 % move all the wave number energy evolution files to that directory
          movefile('Energy_evolutionfig*',mk))                 %movee energy plots to mk
%          movefile('omegafig*',mk)                 %movee omega plots to mk
          
    else
          %movefile('Evolution_subplot*',mk)                 % move all the video files to that directory
          %movefile('filmEvolution*',mk)                 % move all the video files to that directory
          %movefile('Fourier_evolution*',mk)                 % move all the video files to that directory
          %movefile('harmonic_evolution*',mk)                 % move all the image files to that directory
          %movefile('Final_Sk_snapshot*',mk)                 % move all the image files to that directory
          %movefile('Maximum_energy_frequency*',mk)                 % move all the image files to that directory
          %movefile(strcat('*_',num2str(L_flat),'_*.mat'),mk)                 % move the .mat file to that directory
          %movefile('wave_num_fig*',mk)                 % move all the wave number energy evolution files to that directory
          movefile('Energy_evolutionfig*',mk)                 %movee energy plots to mk
          %movefile('omegafig*',mk)                 %movee omega plots to mk
      
    end

end