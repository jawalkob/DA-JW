function [] = rekoTigreReduced(fileName,numAngles, distance, resPath, varargin)
   
    %%Validate inputs
    if ~isstring(fileName)
        error("fileName needs to be of type string")
    end
    
    %if isempty(filePath)
    %    filePathBool = boolean(0);
    %elseif ~isstring(filePath)
    %    error("filePath needs to be of type string")
    %end
    
    if ~isnumeric(numAngles)
        error("numAngles needs to be of type Integer")
    end
    
    if ~isnumeric(distance)
        error("distance needs to be of type numeric")
    end
  
    if isempty(varargin)
        warning("no algorithm choice provided - using FDK")
        algos = ["FDK"];
    else
        algos = string(varargin{:});
      
    end
    
    %%Make Folder for Results
    mkdir(resPath, "RekoRes")
    resPath=resPath+"RekoRes\";

    %% Geometry for Z-position 200 und %170
    zPos = distance;
    switch zPos
        case 200
            % Distances
            geo.DSD = 824.34;                           % Distance Source Detector      (mm)
            geo.DSO = 30.38;                            % Distance Source Origin        (mm)
            % Image parameters
            geo.dVoxel=[0.00468;0.00468;0.00468].*3;                                                   
            % Optional Parameters                                            
            geo.COR=15.520089574216115;                 % y direction displacement for 
                                                        % centre of rotation
                                                        % correction                   (mm)
        case 170
            % Distances
            geo.DSD = 854.34;
            geo.DSO = 60.38;
            % Image parameters
            geo.dVoxel=[0.008976;0.008976;0.008976];
            % Optional Parameters
            geo.COR=8.79141306533552;
    end
    
    % Detector parameters
    geo.nDetector=[768;1066];					        % number of pixels              (px)
    geo.dDetector=[0.127; 0.127].*3; 					    % size of each pixel            (mm)
    geo.sDetector=geo.nDetector.*geo.dDetector;         % total size of the detector    (mm)
    
    % Image parameters
    geo.nVoxel=[1066;1066;768];                        % number of voxels              (vx)
    geo.sVoxel=geo.dVoxel.*geo.nVoxel;
    %geo.sVoxel=[256;256;256];                          % total size of the image       (mm)
    %geo.dVoxel=geo.sVoxel./geo.nVoxel;                 % size of each voxel            (mm)
    
    % Offsets
    geo.offOrigin =[0;0;0];                             % Offset of image from origin   (mm)              
    geo.offDetector=[0; 0];                             % Offset of Detector            (mm)  
                                            
    % Auxiliary 
    geo.accuracy=0.5;                           % Variable to define accuracy of
                                                % 'interpolated' projection
                                                % It defines the amoutn of
                                                % samples per voxel.
                                                % Recommended <=0.5             (vx/sample)
    
    % Optional Parameters
        % There is no need to define these unless you actually need them in your
        % reconstruction
        
        %geo.rotDetector=[0;0;0];                    % Rotation of the detector, by 
                                                        % X,Y and Z axis respectively. (rad)
                                                        % This can also be defined per
                                                         % angle
    %Type
    geo.mode='cone';                            % Or 'parallel'. Geometry type.  
    
    %% Set up angles and load projections:
    
    angles=linspace(0,2*pi,numAngles);
    fileNameBase = append(fileName,"%04d"); % '%04d' bewirkt das in sprintf die Zahl mit Nullen auf 4 Stellen aufgef??llt wird% prealocate
    numFiles = length(angles);

    % Set up workers/pool
    parpool();
    
    % Preallocate    
    proj = zeros(geo.nDetector(1),geo.nDetector(2),numFiles,'single','codistributed');
    proj = gather(proj);
    disp("preallocation done")
    % Load Images    
    parfor j=1:numFiles
        proj(:,:,j)=single(imresize(imread(sprintf(fileNameBase,j-1),"tif"),1/3))
    end
    disp("done loading images")
    % Shut down pool    
    delete(gcp('nocreate'));
    
    %% Invert
    proj=-log(proj/(max(proj(:))+1)); % Beer-Lanbert law
    disp("inversion done")
    %% Reconstruct
    niter=60;

    if any(strcmpi("FDK",algos)) 
        imgFDK=FDK(proj,geo,angles);
        save(resPath+fileName+"IMGFDK" ,"imgFDK",'-v7.3');
        %niftiwrite(imgFDK,append(char(fileName),'FDK.nii'))
        disp("FDK done")
    end
    
    if any(strcmpi("OSSART",algos))
        
        imgOSSART=OS_SART(proj,geo,angles,niter);
        save("C:\Users\jakob\Nextcloud\Diplomarbeit\MATLAB\RekoRes\"+fileName+"OS-SART",'imgOSSART','-v7.3');
        %niftiwrite(imgOSSART,append(char(fileName),'OSSART.nii'))
    end
    if any(strcmpi("CGLS",algos))
        imgCGLS = CGLS(proj,geo,angles,niter);
        save("C:\Users\jakob\Nextcloud\Diplomarbeit\MATLAB\RekoRes\"+fileName+"CGLS","imgCGLS",'-v7.3')
    end
    if any(strcmpi("OS-ASD-POCS",algos))
        epsilon=im3Dnorm(Ax(imgFDK,geo,angles)-proj,'L2')*0.15;
        alpha1=0.002;
        alpha2=0.2;
        imgOSASDPOCS1=OS_ASD_POCS(proj,geo,angles,niter,...
                    'TViter',25,'maxL2err',epsilon,'alpha',alpha1,... % these are very important
                     'Verbose',false,...% less important.
                      'BlockSize',size(angles,2)/10,'OrderStrategy','angularDistance'); %OSC options
        save("C:\Users\jakob\Nextcloud\Diplomarbeit\MATLAB\RekoRes\"+fileName+"OSASDPOCS1","imgOSASDPOCS1","alpha1","epsilon",'-v7.3');
        imgOSASDPOCS2=OS_ASD_POCS(proj,geo,angles,niter,...
                    'TViter',25,'maxL2err',epsilon,'alpha',alpha2,... % these are very important
                     'Verbose',false,...% less important.
                      'BlockSize',size(angles,2)/10,'OrderStrategy','angularDistance'); %OSC options
        save("C:\Users\jakob\Nextcloud\Diplomarbeit\MATLAB\RekoRes\"+fileName+"imgOSASDPOCS2","imgOSASDPOCS2","alpha2","epsilon",'-v7.3');
    end
disp("done")
end
