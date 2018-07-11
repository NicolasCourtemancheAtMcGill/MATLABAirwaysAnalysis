% clear all homemade function
% if not on base computer use clear and clc
clr;
% will be used to count the amount of empty tissues (no MCh_Equilibrate)
brokencounter = 0;
% will be used to count the amount of tissues with MCh_Equilibrate files
stillalive = 0;
% files with less data points (use ctrl+f)
incomplete = 0;
% time counter for testing purposes
tic
% First search into the folders:
% files are all the folders inside the Nicolas dir
% mostly Lung files
% dirFlags creates a logical matrix for the files that are dirs
% subFolders contains the said dirs, using the logical matrix on the whole
% files struct to extract the dirs

files = dir('D:\\Nicolas\');
dirFlags = [files.isdir];
subFolders = files(dirFlags);
% loop going through all the Lung folders in the the Nicolas dir
% folders #1, #2 and #43 are obsolete (not Lung folders)
for Foldernumber = 3:42
    % string for name of the current analysed lung
    Lungname_1 = subFolders(Foldernumber).name;
    % pathname for dir of the Lung's data folder
    nextdir = sprintf('D:\\Nicolas\\%s\\data\\', Lungname_1);
    % Second search into the folders:
    % same process, but delving into the specific tissue files this time
    files_2 = dir(nextdir);
    dirFlags_2 = [files_2.isdir];
    subFolders_2 = files_2(dirFlags_2);
    TracheaCounter = 0;
    MBCounter = 0;
    PBCounter = 0;
    TracheaFile = 0;
    MBFile = 0;
    PBFile = 0;
    
    % the first two folders are once again obsolete
    % i_tnum = index of tissue number
    for i_tnum = 3:length(subFolders_2)
        % string with only the name of the current tissue
        current_tissue = subFolders_2(i_tnum).name;
        % full path for the name of the current tissue with a
        % Mch_Equilibrate_* to find similarly-named files in .mat
        name_withvar = sprintf('D:\\Nicolas\\%s\\data\\%s\\MCh_Equilibrate_*.mat',Lungname_1,current_tissue);
        name_withvar_svc = sprintf('D:\\Nicolas\\%s\\data\\%s\\MCh_shortvcheck*.mat',Lungname_1,current_tissue);
        % 3rd and final search into the folders:
        % compilation of the files with correct naming
        % Files_pre only contains strings with the names of the files, the
        % files themselves are not loaded into this variable
        Files_pre = dir(name_withvar);
        Files_pre_svc = dir(name_withvar_svc);
        % checking if there are actually MCh_Equilibrate type files in the
        % tissue, otherwise going to the next tissue in the lung.
        % If there are none, adds +1 to number of empty tissues,
        % otherwise adds +1 to number of tissues containing a correct file
        %
        % Also checking if the tissue is considered "Bad" by the checker
        clear ispresent;
        ispresent = strfind(current_tissue, ' B');
        if size(ispresent) > 0
            brokencounter = brokencounter + 1;
            continue
        end
        ispresent = strfind(current_tissue, 'Bad');
        
        if size(ispresent) > 0
            brokencounter = brokencounter + 1;
            continue
        end
        ispresent = strfind(current_tissue, 'bad');
        if size(ispresent) > 0
            brokencounter = brokencounter + 1;
            continue
        end
        if size(Files_pre,1) == 0
            brokencounter = brokencounter + 1;
            break
        end
        % Then checking which category the tissue is in. If not any category, will break and kill the tissue.
        TracheaFile = 0;
        MBFile = 0;
        PBFile = 0;
        
        if size(strfind(current_tissue, 'Trachea')) > 0
            TracheaFile = 1;
            TracheaCounter = TracheaCounter + 1;
        end
        if size(strfind(current_tissue, 'MB')) > 0
            MBFile = 1;
            MBCounter = MBCounter + 1;
        end
        if size(strfind(current_tissue, 'PB')) > 0
            PBFile = 1;
            PBCounter = PBCounter + 1;
        end
        if (TracheaFile == 0) && (MBFile == 0) && (PBFile == 0)
            brokencounter = brokencounter + 1;
            continue
        end
        stillalive = stillalive + 1;
        
        
        % If the files have passed the test, assembles a struct with all
        % the MCh files of the tissue in Files_complete.
        % Contains the loaded full files with all the data from the
        % experiments.
        for i_fnum = 1:length(Files_pre)
            Files_complete(i_fnum) = load(sprintf('%s\\%s',Files_pre(i_fnum).folder,Files_pre(i_fnum).name));
            
        end
        %
        % Creating vector with correct order of indexes
        %
        clear i_fnum2;
        for i_fnum2 = 1:length(Files_complete)
            unsortedMoments(i_fnum2) = cellstr(Files_complete(i_fnum2).data.When);
        end
        [sortedMoments, indexes] = sort(unsortedMoments);
        %
        % Removing the equilibrate files that come after the first shortvcheck
        %
        i_stop = -57;
        if size(Files_pre_svc,1) ~= 0
            i_stop2 = length(Files_pre_svc);
            svc_position = 1;
            for i_fnum_svc = 1:i_stop2
                searchforWhen = load(sprintf('%s\\%s',Files_pre_svc(i_fnum_svc).folder,Files_pre_svc(i_fnum_svc).name));
                if isfield(searchforWhen, 'data') ~= 1
                    i_fnum_svc = i_fnum_svc + 1;
                end
                
                if isfield(searchforWhen, 'data') == 1
                    if isfield(searchforWhen.data, 'When') ~= 1
                        i_fnum_svc = i_fnum_svc + 1;
                    end
                    if isfield(searchforWhen.data, 'When') == 1
                        Files_complete_svc(svc_position) = load(sprintf('%s\\%s',Files_pre_svc(i_fnum_svc).folder,Files_pre_svc(i_fnum_svc).name));
                        svc_position = svc_position + 1;
                    end
                end
            end
            for i_fnum2_svc = 1:length(Files_complete_svc)
                unsortedMoments_svc(i_fnum2_svc) = cellstr(Files_complete_svc(i_fnum2_svc).data.When);
            end
            [sortedMoments_svc, ~] = sort(unsortedMoments_svc);
            
            clear i_fnum2
            
            for i_fnum2 = 1:length(Files_complete)
                svcCell(1) = sortedMoments(i_fnum2);
                svcCell(2) = sortedMoments_svc(1);
                [~, svcCell_indexes] = sort(svcCell);
                if svcCell_indexes(1) == 2;
                    i_stop = i_fnum2;
                    clear svcCell svcCell_indexes
                    break
                end
                clear svcCell svcCell_indexes
            end
        end
        if i_stop == -57
            i_stop = length(indexes);
        end
        %
        % Sorting Files_complete into a new vector using the indexes
        %
        clear i_fnum3;
        for i_fnum3 = 1:i_stop
            sorted_Files_complete(i_fnum3) = Files_complete(indexes(i_fnum3)).data;
        end
        %% SMALL VERSION
        %
        %
        from = [1 6000 12000 24000 30000];
        to = [6000 12000 24000 30000 42100];
        smallFrequencies_pre = zeros(length(sorted_Files_complete),5);
        smallAmplitudes_pre = zeros(length(sorted_Files_complete),5);
        
        if length(sorted_Files_complete(1).Fin) < 30001
            
            clear Files_complete Force i sorted_Files_complete unsortedMoments sortedMoments indexes Force2 current_tissue;
            clear name_withvar Files_pre a loc3 v max2 loc2 loc4 peakdistance Frequency means Amplitude max1 path;
            clear name_withvar_svc Files_pre_svc i_fnum_svc Files_complete_svc i_fnum2_svc unsortedMoments_svc sortedMoments_svc;
            clear indexes_svc i_stop ii iii smallFrequencies smallAmplitudes;
            
            if TracheaFile > 0
                TracheaFile = 0;
                TracheaCounter = TracheaCounter - 1;
            end

            if MBFile > 0
                MBFile = 0;
                MBCounter = MBCounter - 1;
            end
            
            if PBFile > 0
                PBFile = 0;
                PBCounter = PBCounter - 1;
            end
           
            incomplete = incomplete + 1;

            continue
        end
        
        for ii = 1: length(sorted_Files_complete)
            for iii = 1:5
                Force01 = sorted_Files_complete(ii).Fin;
                if iii == 5
                    ender = length(Force01);
                else
                    ender = to(iii);
                end
                Force = Force01(from(iii):ender);
                
                [a,~]=xcorr(smooth(Force)-mean(Force));
                % the max(s) of a are the peaks, loc3 contains the locations
                % of the peaks in the Force vector
                [max1, loc3] = max(a);
                searching=true;
                v=loc3;
                % the box
                % sets a "box" of 500 measure units after the first peak
                % searches for the nearest local maximum in that range
                % the distance between the two peaks (peakdistance) is the period
                %
                while searching == true
                    [max2, loc2] = max(a(v:(v+500)));
                    % the peak must be at least 10 units from the edges of the box,
                    % otherwise it cannot be considered a peak
                    if loc2 > 10 & loc2 < 490
                        loc4 = v;
                        searching=false;
                    end
                    v=v+1;
                    if v>length(a)-505
                        searching=false;
                        loc4=loc3;
                    end
                end
                % may be referred to as "period"
                peakdistance = loc2;
                % period of 0 is not a period
                % period above 10k is not a frequency we are looking for
                if peakdistance>0 && peakdistance<10000
                    Frequency = 1/(peakdistance/100);
                    % MOVING AVERAGE
                    % will be used to remove the baseline drift from the Force signal
                    % Force2 is the signal balanced around 0
                    for i=[1:length(Force)-peakdistance]
                        means(i)=mean(Force(i:round(i+peakdistance)));
                        Force2(i)= Force(round(i+0.5*peakdistance))-mean(Force(i:round(i+peakdistance)));
                    end
                    % rms (root mean square) makes the function positive to get the amplitude
                    Amplitude = rms(Force2);
                else
                    Frequency = 0;
                    
                    peakdistance=3500;
                    for i=[1:length(Force)-peakdistance]
                        means(i)=mean(Force(i:round(i+peakdistance)));
                        Force2(i)= Force(round(i+0.5*peakdistance))-mean(Force(i:round(i+peakdistance)));
                    end
                    
                end
                % Second Run for adjusted frequency
                % exact same thing as past block
                % Force2 is re-smoothed to try and get an even better frequency
                % sometimes is useful, sometimes not
                [a,b]=xcorr(smooth(Force2)-mean(Force2));
                [max1, loc3] = max(a);
                searching=true;
                v=loc3;
                peakdistance=0;
                while searching == true
                    
                    [max2, loc2] = max(a(v:(v+500)));
                    if loc2 > 10 & loc2 < 490
                        loc4 = v;
                        searching=false;
                    end
                    v=v+1;
                    if v>length(a)-505
                        searching=false;
                        loc4=loc3;
                    end
                end
                peakdistance = loc2;
                
                if peakdistance > 0
                    Frequency = 1/(peakdistance/100);
                else
                    Frequency = 0;
                    Amplitude = 0;
                end
                Amplitude = rms(Force2);
                
                smallFrequencies_pre(ii, iii) = Frequency;
                smallAmplitudes_pre(ii, iii) = Amplitude;
                clear Force01 Force a max1 loc3 v max2 loc2 loc4 peakdistance Frequency means i Force2 Amplitude
            end
        end
        rows = size(smallFrequencies_pre,1);
        if size(smallFrequencies_pre, 1) ~= 1
            smallFrequencies = sum(smallFrequencies_pre)/rows;
            smallAmplitudes = sum(smallAmplitudes_pre)/rows;
        else
            smallFrequencies = smallFrequencies_pre;
            smallAmplitudes = smallAmplitudes_pre;
        end
        
        
        %% COMPLETE VERSION
        %
        % If there is more than one MCh file, appending them together to
        % make the "Force" file. Otherwise Force will only contain one.
        %
        if length(sorted_Files_complete) > 1
            Force=sorted_Files_complete(1).Fin;
            for i = 2:length(sorted_Files_complete)
                Force = [Force sorted_Files_complete(i).Fin];
            end
        end
        if length(sorted_Files_complete) == 1
            Force = sorted_Files_complete(1).Fin;
        end
        %% FORCE -> FORCE2
        %
        % Removing the mean of the Force
        % smooth helps but not necessary
        % xcorrelating to find where the peaks are
        %
        [a,~]=xcorr(smooth(Force)-mean(Force));
        % the max(s) of a are the peaks, loc3 contains the locations
        % of the peaks in the Force vector
        [max1, loc3] = max(a);
        searching=true;
        v=loc3;
        % the box
        % sets a "box" of 500 measure units after the first peak
        % searches for the nearest local maximum in that range
        % the distance between the two peaks (peakdistance) is the period
        %
        while searching == true
            [max2, loc2] = max(a(v:(v+500)));
            % the peak must be at least 10 units from the edges of the box,
            % otherwise it cannot be considered a peak
            if loc2 > 10 & loc2 < 490
                loc4 = v;
                searching=false;
            end
            v=v+1;
            if v>length(a)-505
                searching=false;
                loc4=loc3;
            end
        end
        % may be referred to as "period"
        peakdistance = loc2;
        % period of 0 is not a period :)
        % period above 10k is not a frequency we are looking for
        if peakdistance > 0 & peakdistance < 10000
            Frequency = 1/(peakdistance/100);
            % MOVING AVERAGE
            % will be used to remove the baseline drift from the Force signal
            % Force2 is the signal balanced around 0
            for i=[1:length(Force)-peakdistance]
                means(i)=mean(Force(i:round(i+peakdistance)));
                Force2(i)= Force(round(i+0.5*peakdistance))-mean(Force(i:round(i+peakdistance)));
            end
            % rms (root mean square) makes the function positive to get the amplitude
            Amplitude = rms(Force2);
        else
            Frequency = 0;
            
            peakdistance=3500;
            for i=[1:length(Force)-peakdistance]
                means(i)=mean(Force(i:round(i+peakdistance)));
                Force2(i)= Force(round(i+0.5*peakdistance))-mean(Force(i:round(i+peakdistance)));
            end
            
        end
        % Second Run for adjusted frequency
        % exact same thing as past block
        % Force2 is re-smoothed to try and get an even better frequency
        % sometimes is useful, sometimes not
        [a,b]=xcorr(smooth(Force2)-mean(Force2));
        [max1, loc3] = max(a);
        searching=true;
        v=loc3;
        peakdistance=0;
        while searching == true
            
            [max2, loc2] = max(a(v:(v+500)));
            if loc2 > 10 & loc2 < 490
                loc4 = v;
                searching=false;
            end
            v=v+1;
            if v>length(a)-505
                searching=false;
                loc4=loc3;
            end
        end
        peakdistance = loc2;
        
        if peakdistance>0
            Frequency = 1/(peakdistance/100);
        else
            Frequency = 0;
            Amplitude = 0;
        end
        Amplitude = rms(Force2);
        %% FFT VARIABLES
        % This section is used to create the necessary variables to plot an FFT
        % This can be used to check (afterwards not), if the function is paused,
        % if the frequency that was found exists in the spectrum
        Fs = 100;            % samples per second/sampling frequency
        L = length(Force2); % Length of signal
        Y = fft(Force2);        % transform of input signal
        %Code from matlab website on fourier transforms
        P2 = abs(Y/L);                %
        P1 = P2(1:L/2+1);             %
        P1(2:end-1) = 2*P1(2:end-1);  %
        f = Fs*(0:(L/2))/L;           %
        % %%
        % Plot section (to do manually)
        % plot(f,P1)                                        % |Until the end: plotting information
        % title('Single-Sided Amplitude Spectrum of X(t)')  % |
        % xlabel('f (Hz)')                                  % |
        % ylabel('|P1(f)|')                                 % |
        
        %% PLOTTING
        %
        % correct variables to plot the graph of the final force signal.
        %
        path = sprintf('%s',Files_pre(i_fnum).folder);
        h=figure(1);
        plot([1:length(Force2)]/6000,Force2);
        xlabel('Time(min)');
        ylabel('Force(mN)');
        line([0 length(Force2)/6000],[Amplitude Amplitude],'color',[0 0 0]);
        if Frequency==0
            text(0.1*length(Force2)/6000,1.2*max(Force2),['Frequency = ' num2str(Frequency)])
            text(0.1*length(Force2)/6000,0.8*max(Force2),['RMS Amplitude = ' num2str(Amplitude)])
        else
            text(0.1*length(Force2)/6000,2*Amplitude,['Frequency = ' num2str(Frequency)])
            text(0.1*length(Force2)/6000,1.3*Amplitude,['RMS Amplitude = ' num2str(Amplitude)])
        end
        Dummy=strsplit(path,'\\');
        Dummy2=strsplit(Dummy{3},' ');
        Lungno=[(Dummy2(1:2))];
        TitleArray={'Drift Corrected Force',[Lungno{1} ' ' Lungno{2} ' ' Dummy{5}]};
        title(TitleArray);
        %
        % saving the graph
        %
        saveas(h,['C:\\Nicolas\\LungDataTesting_4\\' TitleArray{2} '.png']);
        
        if TracheaFile > 0
            FinalFrequenciesTrachea_pre(TracheaCounter, 1:5) = smallFrequencies(1:5);
            FinalAmplitudesTrachea_pre(TracheaCounter, 1:5) = smallAmplitudes(1:5);
        end
        if MBFile > 0
            FinalFrequenciesMB_pre(MBCounter, 1:5) = smallFrequencies(1:5);
            FinalAmplitudesMB_pre(MBCounter, 1:5) = smallAmplitudes(1:5);
        end
        if PBFile > 0
            FinalFrequenciesPB_pre(PBCounter, 1:5) = smallFrequencies(1:5);
            FinalAmplitudesPB_pre(PBCounter, 1:5) = smallAmplitudes(1:5);
        end
        %
        % Clearing for next loop
        % THIS SHOULD BE AT THE END OF THE LOOP
        clear Files_complete Force i sorted_Files_complete unsortedMoments sortedMoments indexes Force2 current_tissue;
        clear name_withvar Files_pre a loc3 v max2 loc2 loc4 peakdistance Frequency means Amplitude max1 path;
        clear name_withvar_svc Files_pre_svc i_fnum_svc Files_complete_svc i_fnum2_svc unsortedMoments_svc sortedMoments_svc;
        clear indexes_svc i_stop ii iii smallFrequencies smallAmplitudes;
    end
    % //
    % //end of the Tissue loop//
    % //
    % SAVING
    % Trachea
    if TracheaCounter > 0
        if size(FinalFrequenciesTrachea_pre, 1) ~= 1
            FinalFrequenciesTrachea = sum(FinalFrequenciesTrachea_pre)/TracheaCounter;
            FinalAmplitudesTrachea = sum(FinalAmplitudesTrachea_pre)/TracheaCounter;
        else
            FinalFrequenciesTrachea = FinalFrequenciesTrachea_pre;
            FinalAmplitudesTrachea = FinalAmplitudesTrachea_pre;
        end
        [~,~,C]=xlsread('C:\\Nicolas\\DataTesting_4.xlsx');
        %         -----------freq storing
        C(end+1, 1:7) = {[Dummy{3} ' ' Dummy{5}], 'Trachea', FinalFrequenciesTrachea(1), FinalFrequenciesTrachea(2), FinalFrequenciesTrachea(3), FinalFrequenciesTrachea(4), FinalFrequenciesTrachea(5)};
        %         -----------amp storing
        C(end, 8:12) = {FinalAmplitudesTrachea(1), FinalAmplitudesTrachea(2), FinalAmplitudesTrachea(3), FinalAmplitudesTrachea(4), FinalAmplitudesTrachea(5)};

        xlswrite('C:\\Nicolas\\DataTesting_4.xlsx',C);
    end
    % MB
    if MBCounter > 0
        if size(FinalFrequenciesMB_pre, 1) ~= 1
            FinalFrequenciesMB = sum(FinalFrequenciesMB_pre)/MBCounter;
            FinalAmplitudesMB = sum(FinalAmplitudesMB_pre)/MBCounter;
        else
            FinalFrequenciesMB = FinalFrequenciesMB_pre;
            FinalAmplitudesMB = FinalAmplitudesMB_pre;
        end
        [~,~,C]=xlsread('C:\\Nicolas\\DataTesting_4.xlsx');
        %         -----------freq storing
        C(end+1, 1:7) = {[Dummy{3} ' ' Dummy{5}], 'MB', FinalFrequenciesMB(1), FinalFrequenciesMB(2), FinalFrequenciesMB(3), FinalFrequenciesMB(4), FinalFrequenciesMB(5)};
        %         -----------amp storing
        C(end, 8:12) = {FinalAmplitudesMB(1), FinalAmplitudesMB(2), FinalAmplitudesMB(3), FinalAmplitudesMB(4), FinalAmplitudesMB(5)};

        xlswrite('C:\\Nicolas\\DataTesting_4.xlsx',C);
    end
    % PB
    if PBCounter > 0
        if size(FinalFrequenciesPB_pre, 1) ~= 1
            FinalFrequenciesPB = sum(FinalFrequenciesPB_pre)/PBCounter;
            FinalAmplitudesPB = sum(FinalAmplitudesPB_pre)/PBCounter;
        else
            FinalFrequenciesPB = FinalFrequenciesPB_pre;
            FinalAmplitudesPB = FinalAmplitudesPB_pre;
        end
        [~,~,C]=xlsread('C:\\Nicolas\\DataTesting_4.xlsx');
        %         -----------freq storing
        C(end+1, 1:7) = {[Dummy{3} ' ' Dummy{5}], 'PB', FinalFrequenciesPB(1), FinalFrequenciesPB(2), FinalFrequenciesPB(3), FinalFrequenciesPB(4), FinalFrequenciesPB(5)};
        %         -----------amp storing
        C(end, 8:12) = {FinalAmplitudesPB(1), FinalAmplitudesPB(2), FinalAmplitudesPB(3), FinalAmplitudesPB(4), FinalAmplitudesPB(5)};

        xlswrite('C:\\Nicolas\\DataTesting_4.xlsx',C);
    end
    clear TracheaCounter TracheaFile FinalAmplitudesTrachea FinalFrequenciesTrachea FinalFrequenciesTrachea_pre FinalAmplitudesTrachea_pre
    clear MBCounter MBFile FinalAmplitudesMB FinalFrequenciesMB FinalFrequenciesMB_pre FinalAmplitudesMB_pre
    clear PBCounter PBFile FinalAmplitudesPB FinalFrequenciesPB FinalFrequenciesPB_pre FinalAmplitudesPB_pre
    clear C
end
% //
% //end of the Lung loop//
% //

toc
fprintf('\nNumber of tissues containing no MCh_Equilibrate files: %i\n', brokencounter);
fprintf('\nNumber of tissues containing MCh Equilibrate files: %i\n', stillalive);
