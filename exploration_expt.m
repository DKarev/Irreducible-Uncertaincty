function EXPT = exploration_expt(local)

    % creates EXPT structure for CCNL fMRI processing pipeline
    %
    % USAGE: EXPT = contextExpt()
    %
    % INPUTS:
    %   local (optional) - true if file path is on local computer, false if on NCF
    %
    % OUTPUTS:
    %   EXPT - experiment structure with fields
    %          .TR - repetition time
    %          .create_multi - function handle for creating multi structure
    %          .modeldir - where to put model results
    %          .subject(i).datadir - directory for subject data
    %          .subject(i).functional - .nii files for runs
    %          .subject(i).structural - .nii for structural scan
    %
    % Cody Kommers, July 2016
    % Momchil Tomov, Nov 2016
    
    % set default parameters
    %
    if nargin < 1
        [~, name] = system('hostname');
        if ~isempty(strfind(name, 'omchil'))
            % err on the side of falsely thinking it's NCF. Because locally
            % you will catch that mistake immediatley. On NCF, you will
            % catch it after you're already sent 100 jobs and they all
            % fail 2 days later...
            %
            local = true;
        else
            local = false;
        end
    end
    
    % set main directory
    %
    if local
        %exptdir = '/Users/momchil/Dropbox/research/context/'; % local group level on dropbox
       % exptdir = '/Volumes/fMRI/Exploration/'; % local group level on external drive
        exptdir = '/Users/momchil/Dropbox/Research/exploration/';
    else
        exptdir = '/ncf/gershman/Lab/Exploration/'; % on CBS central server
    end
    
   % % Load data from file with all subjects
   % %
   % [data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
   % 
    [allSubjects, subjdirs, goodRuns] = exploration_getSubjectsDirsAndRuns();
   % allSubjects
   % metadata.allSubjects
   % assert(isequal(allSubjects, metadata.allSubjects));
    
    for subj = 1:length(allSubjects)
        subjdir = [exptdir, 'subjects/', subjdirs{subj}, '/'];
        EXPT.subject(subj).datadir = [subjdir, 'preproc'];
        
        EXPT.subject(subj).structural = 'struct.nii';
        
        
        %assert(nRuns{subj} == length(unique(data.runId(strcmp(data.participant, allSubjects{subj})))));
        
       
        i = 1;
        for run = 1:length(goodRuns{subj})
            if(goodRuns{subj}(run))
                EXPT.subject(subj).functional{i} = ['run',sprintf('%03d',run),'.nii'];
                i = i + 1;
            end 
        end
        disp(EXPT.subject(subj));
    end
    
    % TR repetition time
    EXPT.TR = 2; %seconds

    % Function handle to create subject multi structure
    EXPT.create_multi = @exploration_create_multi;
    % Function handle to create subject RSA structure
    EXPT.create_rsa = @exploration_create_rsa;

    % Where you want model output da:ta to live
    EXPT.modeldir = [exptdir, 'glmOutput'];
    % Where the RSA output lives
    EXPT.rsadir = [exptdir, 'rsaOutput'];
    
    % Where the data live, but not sure which data
    EXPT.datadir = [exptdir, 'testOutput'];


end
