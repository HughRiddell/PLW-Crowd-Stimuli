function [] = walker_field_dots_final()

% experiment to measure heading performance in crowds of normal PLWs that
% articulate their limbs or maintain a static posture while translating
% through space. Ground plane points are also presented and varied.
% Written by Hugh Riddell, 2017

%-------------
% Input initial variables
%-------------

ID = input('Enter subject ID ','s'); %Input subject ID, this will also be the file name of the output
practice = input('Practice run [1] Experimental run [0] '); %input whether this is a practice run [1] or not [0]
static = input('Static [0] or articulating walkers [1]? '); % whether or not walkers maintain a static pose [0] or articulate [1] while translating

%-------------
% Set up stuff
%-------------

% GL data structure needed for all OpenGL demos:
global GL;

% Is the script running in OpenGL Psychtoolbox? Abort, if not.
AssertOpenGL;

% Restrict KbCheck to checking of ESCAPE key:
KbName('UnifyKeynames'); 
exitkey=KbName('ESCAPE');

%skip screen tests stops Psychtoolbox from crashing
Screen('Preference','Verbosity',1);
Screen('Preference', 'SkipSyncTests', 1);

% Find the screen to use for display:
screenid=max(Screen('Screens'));
stereoMode = 0;
multiSample = 0;

% get resolution in px
scr=Screen('Resolution',screenid);
pix_diag=sqrt(scr.width^2+scr.height^2);

%physical screen dimensions in cm. 
%Note: must be adjusted for each individual display if not set for the
%correct display stimuli will not be correct!

%Currently set up to work on Macs using display at Muenster. Some unresolved issues
%when porting to Windows laptop (Jan 2018).

scr.pwidth=248; %screen physical width in cm
scr.pheight=182; %screen physical height in cm
scr.diag = sqrt(scr.pwidth^2+scr.pheight^2); %screen physical size on diagonal

pix_per_cm = pix_diag/scr.diag; %number of px per cm

distance1=100; %physical viewing distance in cm


grad = @(x) (tand(x/2)*(distance1*2))*pix_per_cm; % inline function to calculate pixel position of objects. input the degrees, fucntion outputs px conversion

%-----------
% Parameters for experiment
%-----------

fps = 60; % frame rate of monitor, input manually. PTB fuction can output funny values on LCD monitors

nframes = 120;  %duration of stimulus, in frames. MS display operates at 60fps.

ntrials = 15;    %number of repartitions for each condition

numwalkers = 8; %number of walkers in crowd

d=10;   %scene depth in meters

hdrange = 12; %heading range in degrees

tspeed=1.5/fps;  %speed which the observer translates through the environment 

if practice 
    
    ntrials=15; %no. practice trials 30
    
end

conditions = [0, ceil(exp(linspace(log(1),log(20),8)))]; % condition with different numbers of fludicial (ground) points. pseudo log spaced between 0 and 20

conditions=Shuffle(conditions); % randomise conditions for each block

%-----------
% OpenGL Setup
%-----------

% Setup Psychtoolbox for OpenGL 3D rendering support and initialize the
% mogl OpenGL for Matlab wrapper:
InitializeMatlabOpenGL;

PsychImaging('PrepareConfiguration');

% Open a double-buffered full-screen window on the main displays screen.
[win, winRect] = PsychImaging('OpenWindow', screenid, 0, [], [], [], stereoMode, multiSample);
[xaxis, yaxis] = RectCenter(winRect); %center of window

HideCursor; % hide mouse
Priority(MaxPriority(win)); %set highest computing priority

% Setup the OpenGL rendering context of the onscreen window for use by
% OpenGL wrapper. After this command, all following OpenGL commands will
% draw into the onscreen window 'win':
Screen('BeginOpenGL', win);

% Set viewport properly:
glViewport(0, 0, scr.width, scr.height);

% Setup OpenGL local lighting model: The lighting model supported by
% OpenGL is a local Phong model with Gouraud shading.

% Enable the first local light source GL.LIGHT_0. Each OpenGL
% implementation is guaranteed to support at least 8 light sources,
% GL.LIGHT0, ..., GL.LIGHT7
glEnable(GL.LIGHT0);

% Enable alpha-blending for smooth dot drawing:
glEnable(GL.BLEND);
glBlendFunc(GL.SRC_ALPHA, GL.ONE_MINUS_SRC_ALPHA);

glEnable(GL.DEPTH_TEST);

% Set projection matrix: This defines a perspective projection,
% corresponding to the model of a pin-hole camera - which is a good
% approximation of the human eye and of standard real world cameras --
% well, the best aproximation one can do with 3 lines of code ;-)
glMatrixMode(GL.PROJECTION);
glLoadIdentity;

% set up perspective projection
gluPerspective(89, scr.width/scr.height, 0.5, d);

% Setup modelview matrix: This defines the position, orientation and
% looking direction of the virtual camera:
glMatrixMode(GL.MODELVIEW);
glLoadIdentity;

% Our point lightsource is at position (x,y,z) == (1,2,3)...
glLightfv(GL.LIGHT0,GL.POSITION,[ 1 2 3 0 ]);

% Set background clear color to 'black' (R,G,B,A)=(0,0,0,0):
glClearColor(0,0,0,0);

% Clear out the backbuffer: This also cleans the depth-buffer for
% proper occlusion handling: You need to glClear the depth buffer whenever
% you redraw your scene, e.g., in an animation loop. Otherwise occlusion
% handling will screw up in funny ways...
glClear(GL.DEPTH_BUFFER_BIT);

% Finish OpenGL rendering into PTB window. This will switch back to the
% standard 2D drawing functions of Screen and will check for OpenGL errors.

Screen('EndOpenGL', win);

% Show rendered image at next vertical retrace:
Screen('Flip', win);

%-----------
% start of experiment
%-----------

Screen('TextSize',win, 20);

experiment_trials = 1; %initialize trial counter   

%loop for each condition block
for select_cond = 1:length(conditions)
    
    ndots = conditions(1,select_cond); %select number of dots for cond.

    %display welcome message
    [~, ~, buttons1]=GetMouse(screenid);
    Screen('TextSize',win, 20);
    white = WhiteIndex(win);
    
    while ~any(buttons1)
        
        Screen('DrawText',win, 'Click the mouse to begin the block.',xaxis-200,yaxis,white);
        Screen('DrawingFinished', win);
        Screen('Flip', win);
        [~, ~, buttons1]=GetMouse(screenid);
        
    end

    %set/reset observer position in scene
    translate=0; 

    %save origin directory where files are
    origin = pwd;
    
    FID = fopen('sample_walker3.txt');    %open walker data file
     
    walker_array = fscanf(FID,'%f');     %read into matlab
    fclose(FID);

    walker_array=reshape(walker_array,3,[]).*0.00001;  %order and scale walker array   

    
    %-----------
    % Trial Loop
    %-----------
    
    %loop through ntrials for each condition. Conditions are blocked by
    %number of ground points
    
    for trials = 1:ntrials     
  
    
    %-----------
    % set up walkers
    %-----------   
    
    %randomly select starting phase
    clear walker_phase    
    numorder=(1:16:length(walker_array));
    walker_phase(1:numwalkers)=datasample(numorder,numwalkers,'Replace',true);
     
    walk=zeros(1,numwalkers); %set walked distance to zero for all walkers
    walk_speed = 0.01; % walking speed for walkers

    %set walker rotation for each walker
    walker_rotation = 360*rand(1,numwalkers);

    %build rotation matrix for walkers this provides a way to randomly orient
    %each walker in 3D
    for build_r_mat = 1:numwalkers
        
        r{build_r_mat} = [cosd(walker_rotation(build_r_mat)), 0, -sind(walker_rotation(build_r_mat));...
                         0,                                   1, 0;...
                         sind(walker_rotation(build_r_mat)),  0, cosd(walker_rotation(build_r_mat))];
                     
    end

    %gen walker random starting positions
    [walkerX,walkerY,walkerZ] = CreateUniformDotsIn3DFrustum(numwalkers,89,scr.width/scr.height,0.5,d,1.4); %generate walker positions randomly within view frustum
    walkerX = linspace(-3,3,numwalkers)+rand(1,numwalkers).*0.5.*datasample([1 -1],numwalkers,'Replace',true); %restrict x positions to  make a denser crowd

     
    %-----------
    % set up heading
    %----------- 
    
    Screen('BeginOpenGL',win)
    
    glLoadIdentity
    
    viewport=glGetIntegerv(GL.VIEWPORT); %get viewport
    modelview=glGetDoublev(GL.MODELVIEW_MATRIX); %get modelview matrix            
    projection=glGetDoublev(GL.PROJECTION_MATRIX); %get projection matrix

    head = hdrange*rand(); %get random heading position 
    randheadpos=grad(head); %convert to px
    randhead = datasample([1 -1],1,'Replace',true); %add random L(+)/R(-) position so that heading is left or right of center

    [heading_world, trash, trash]=gluUnProject(xaxis+randheadpos*randhead,0,1,modelview,projection,viewport); %get position of heading in OpenGL coords
    
    %gen random dot positions for ground plane
    [dotx, doty, dotz] = CreateUniformDotsIn3DFrustum(ndots,89,scr.width/scr.height,0.5,d,1.4);
    
    % if one dot condition restrict position so dot is not in periphery
    if ndots == 1
        
        dotx = (3+3).*rand()-3;
        doty = -1.4;
        dotz = ((d-d*0.5).*rand(1,ndots)+(0.5*d)).*-1;
        
    end
    
    dot_xyz = [dotx;doty;dotz]; %create xyz-matrix for ground plane dots  
    
    %-----------
    % set up view frustum culling
    %-----------   
    
    %frustum culling is used later to identify walkers/dots that have left
    %FOV
    
    vangle=atand((heading_world)/(d)); %viewing angle of simulated observer
     
    %rotation matrix tp get sim observer looking direction
    rdots = [cosd(-vangle), 0, -sind(-vangle);...
            0,             1, 0;...
            sind(-vangle),  0, cosd(-vangle)];
        
    glPushMatrix   
    
        glLoadIdentity

        glRotatef(-vangle,0,1,0) %rotate to emulate viewing angle

        proj=glGetFloatv(GL.PROJECTION_MATRIX); %get proj matrix
        modl=glGetFloatv(GL.MODELVIEW_MATRIX); %get modelview matrix
   
    glPopMatrix
    
    %make matrices 4x4
    modl=reshape(modl,4,4); 
    proj=reshape(proj,4,4);

    frustum=getFrustum(proj,modl); %function to get view frustum

    Screen('EndOpenGL', win)
   
    %-----------
    % Animation Loop
    %-----------   

    % draws each frame for crowd animation
    
    for i = 1:nframes;

        Screen('BeginOpenGL',win);

        glLoadIdentity

        gluLookAt(0,0,0,heading_world,0,-d,0,1,0); %set camera to look at heading while located at center, heading angle is defined by the difference between these two  

        glTranslatef(0,0,translate) %translate scene - observer motion

        %-----------
        % Fludicial point drawing
        %-----------   

        %culling and repositioning for dots that have left the view
        
        for test = 1:ndots

            %get point for culling
            p1 = [dot_xyz(1,test),dot_xyz(2,test),dot_xyz(3,test)+translate]*rdots; 

            %normalize
            p1=p1/norm(p1);

            %test and cull
            if  frustum(1,1)*p1(1) + frustum(1,2)*p1(2) + frustum(1,3)*p1(3) + frustum(1,4) < 0 || frustum(5,1)*p1(1) + frustum(5,2)*p1(2) + frustum(5,3)*p1(3) + frustum(5,4) < 0 || frustum(2,1)*p1(1) + frustum(2,2)*p1(2) + frustum(2,3)*p1(3) + frustum(2,4) < 0 

                [dot_xyz(3,test)] = -d-translate; %reset dots at back of frusutm

            end     
        end

        %these variables set up some point drawing
        nrflowdots=size(dot_xyz,2);
        nvcflow=size(dot_xyz,1);

        glClear(GL.DEPTH_BUFFER_BIT)

        %this bit of code was taken out of the moglDrawDots3D psychtoolbox function which is EXTREMELY inefficient. it is much quicker to just use the relevant openGL function to draw points
        glVertexPointer(nvcflow, GL.DOUBLE, 0, dot_xyz);
        glEnableClientState(GL.VERTEX_ARRAY);

        glEnable(GL.POINT_SMOOTH); %enable anti-aliasing
        glHint(GL.POINT_SMOOTH_HINT, GL.DONT_CARE); %but it doesnt need to be that fancy. they are just white dots after all

        glPushMatrix

            glRotatef(-vangle,0,1,0) 
            glColor3f(0.6,0.6,0.6)
            glPointSize(4)

            glDrawArrays(GL.POINTS, 0, nrflowdots); %draw the fludicial points

        glPopMatrix

        %-----------
        % Draw Walkers
        %-----------      
        
        %loops through and draws each individual walker. Loop is necessary
        %because each walker has unique position, rotation and phase.

        for walker = 1:numwalkers 
            
            %reset phase if we are at the end of the file
            if walker_phase(walker)+16+12 > length(walker_array) 
                
                walker_phase(walker)=1;
                
            end

            %get walker array for frame
            xyzmatrix = walker_array(:,walker_phase(walker):walker_phase(walker)+11);

            %advance frame if the walker is articulating
            if static
                
                walker_phase(:,walker) = walker_phase(:,walker) + 16;
                
            end

            %-----------
            % Reposition lost walker using frustum culling
            %-----------   

            %get walker position
            p = [walkerX(walker)+walk(walker),walkerY(walker),walkerZ(walker)+translate]+[walk(walker),0,0]*r{walker}; 

            %normalize
            p=p/norm(p);

            %test and cull
            if  frustum(1,1)*p(1) + frustum(1,2)*p(2) + frustum(1,3)*p(3) + frustum(1,4) < 0 || frustum(5,1)*p(1) + frustum(5,2)*p(2) + frustum(5,3)*p(3) + frustum(5,4) < 0 || frustum(2,1)*p(1) + frustum(2,2)*p(2) + frustum(2,3)*p(3) + frustum(2,4) < 0 

                walkerZ(walker)= -d-translate; %compensate for moving in depth

                walk(walker)= 0; %reset walked disatance

            end      

            %-----------
            % PLW drawing
            %-----------   

            %these variables set up some point drawing same as for
            %fludicial pts
            
            nrdots=size(xyzmatrix,2);
            nvc=size(xyzmatrix,1);

            glClear(GL.DEPTH_BUFFER_BIT)

            glVertexPointer(nvc, GL.DOUBLE, 0, xyzmatrix);
            glEnableClientState(GL.VERTEX_ARRAY);

            glEnable(GL.POINT_SMOOTH); 
            glHint(GL.POINT_SMOOTH_HINT, GL.DONT_CARE); 
            
            glPushMatrix
            
                glTranslatef(walkerX(walker),walkerY(walker),walkerZ(walker)); %move the points to the right location  

                %do rotation and walking translation
                glRotatef(walker_rotation(walker),0,1,0);
                glTranslatef(walk(walker),0,0);
                walk(walker)=walk(walker)+walk_speed; %update walked position

                glColor3f(0.6,0.6,0.6)
                glPointSize(4)
                glDrawArrays(GL.POINTS, 0, nrdots); %draw the points

            glPopMatrix

        end

        Screen('EndOpenGL',win);

        translate=translate+tspeed; % update observer position

        Screen('Flip', win);

    end    %end animation loop
      
    %-----------
    % Post-animation
    %-----------   
       
    Screen('BeginOpenGL',win)
    
    %reset opengl world
    glMatrixMode(GL.MODELVIEW)
    glLoadIdentity

    Screen('EndOpenGL',win)

    SetMouse(rand()*scr.width,yaxis,win);   %set the mouse in  random location on the screen

    buttons = 0;
        
    %-----------
    % Subject response loop
    %-----------  
    
    while ~buttons

        %abort program early if desired
        [~,~,keycode]=KbCheck();

        if keycode(exitkey)   

            Screen('CloseAll');
            return  

        end

        [mx, ~, buttons]=GetMouse(screenid); %get mouse position on screen

        Screen('DrawLine', win, [255 0 0 0], mx, yaxis+50, mx, yaxis); % draws estimating line on horizon of ground plane
        Screen('DrawingFinished',win);


        %-----------
        % Redraw Scene
        %-----------   

        Screen('BeginOpenGL',win);

        glMatrixMode(GL.MODELVIEW)
        glLoadIdentity
        glClear(GL.DEPTH_BUFFER_BIT)

        %set camera looking position and location
        gluLookAt(0,0,0,heading_world,0,-d,0,1,0);

        glTranslatef(0,0,translate-tspeed);

        %dot drawing again

        nrflowdots=size(dot_xyz,2);
        nvcflow=size(dot_xyz,1);

        glClear(GL.DEPTH_BUFFER_BIT)

        glVertexPointer(nvcflow, GL.DOUBLE, 0, dot_xyz);
        glEnableClientState(GL.VERTEX_ARRAY);

        glEnable(GL.POINT_SMOOTH); 
        glHint(GL.POINT_SMOOTH_HINT, GL.DONT_CARE); 
        
        glColor3f(0.6,0.6,0.6)
        glPointSize(4)

        glPushMatrix
        glRotatef(-vangle,0,1,0)
        glDrawArrays(GL.POINTS, 0, nrflowdots); %draw the points
        glPopMatrix

        % walker drawing again
        for walker = 1:numwalkers

            xyzmatrix = walker_array(:,walker_phase(walker):walker_phase(walker)+11);

            nrdots=size(xyzmatrix,2);
            nvc=size(xyzmatrix,1);

            glClear(GL.DEPTH_BUFFER_BIT)
            
            glVertexPointer(nvc, GL.DOUBLE, 0, xyzmatrix);
            glEnableClientState(GL.VERTEX_ARRAY);

            glEnable(GL.POINT_SMOOTH); 
            glHint(GL.POINT_SMOOTH_HINT, GL.DONT_CARE); 
            
            glPushMatrix
            glTranslatef(walkerX(walker),walkerY(walker),walkerZ(walker)); 
            
           
            glRotatef(walker_rotation(walker),0,1,0);
            glTranslatef(walk(walker),0,0);                    


            glColor3f(0.6,0.6,0.6)
            glPointSize(4)
            glDrawArrays(GL.POINTS, 0, nrdots); %draw points
            glPopMatrix


        end


        Screen('EndOpenGL',win);

        %-----------
        % Calculate heading error
        %-----------   
        
        %stop on mouse click
        if any(buttons)

            %mouse
            
            [mousex, ~, ~] = GetMouse(screenid); %get mouse position in px
            
            mouse_cm = mousex/pix_per_cm; %mouse in cm
            
            va_mouse = 2*atand(mouse_cm/(distance1*2)); %mouse in va
            
            va_mouse_centred = 2*atand(((mousex-xaxis)/pix_per_cm)/(distance1*2)); %mouse in va centered with respect to the screen

            %heading 
            
            headx = xaxis-randheadpos*randhead; %heading positin in px
            
            head_cm = headx/pix_per_cm; %heading in cm 
            
            va_head = 2*atand(head_cm/(distance1*2));  %heading in visual angle
            
            va_head_centred = 2*atand(((headx-xaxis)/pix_per_cm)/(distance1*2)); %heading in visual angle centered with respect to the screen

            %heading error (difference between estimated heading and actual heading positions)
            error_va=va_head-va_mouse;  %in VA
            error_cm=head_cm-mouse_cm;  %in CM

        end

        Screen('Flip', win);

    end

    Screen('Flip',win);

    WaitSecs(0.5); %inter stimulus interval 500ms

    translate = 0; %reset translation counter

    %-----------
    % update and save output variable
    %-----------  

    error_var(experiment_trials,1) = mousex; %mouse pos in px
    error_var(experiment_trials,2) = headx; %heading in px
    error_var(experiment_trials,3) = mouse_cm;  %mouse pos in cm
    error_var(experiment_trials,4) = head_cm; %heading pos in cm
    error_var(experiment_trials,5) = va_mouse;  %mouse pos in va
    error_var(experiment_trials,6) = va_head;   %heading pos in va
    error_var(experiment_trials,7) = va_mouse_centred;  %mouse pos in VA centered with respect to screen
    error_var(experiment_trials,8) = va_head_centred;  %heading pos in VA centered with respect to screen    
    error_var(experiment_trials,9) = error_va; %heading error in VA signed to denote direction (+) = left of center (-) = right of center
    error_var(experiment_trials,10) = error_cm; %heading error in cm signed to denote direction (+) = left of center (-) = right of center
    error_var(experiment_trials,11) = head; % inputed heading position 
    error_var(experiment_trials,12) = randheadpos*randhead; %inputed heading position plus L(+) R(-) direction from center
    error_var(experiment_trials,13) = ndots; %number of static points
    error_var(experiment_trials,14) = static; %whether or not walkers were articulating(1) or static(0)

    experiment_trials = experiment_trials + 1; % update trial counter
        
    %if not practice write a temporary file with current data
    if ~practice

        walker_label = 'norm';

        switch static
            
            case 1
                stat_label = '_articulating';
                
            case 0
                stat_label = '_static';
                
        end

        %save data
        cd('data');

        out = num2cell(error_var);
        subname{experiment_trials-1,1} = ID; %add ID to output
        final = cat(2,out,subname);
        fname=[ID, '_', walker_label, stat_label, '_translating_walker_crowd_repeat_experiment_dots.txt'];
        formatspec='%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%s\n';
        fid = fopen(fname, 'w') ;

        for row = 1:size(final,1);
            fprintf(fid,formatspec,final{row,:});
        end

        fclose(fid);

        cd(origin)

    end

    end
    
    %-----------
    % End Program and save final data
    %-----------  
    
    if ~practice

        walker_label = 'norm';

        switch static
            
            case 1
                stat_label = '_articulating';
                
            case 0
                stat_label = '_static';
                
        end

        %save data
        cd('data');

        out = num2cell(error_var);
        subname{experiment_trials-1,1} = ID;
        final = cat(2,out,subname);
        formatspec='%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%s\n';
        fid = fopen(fname, 'w') ;

        for row = 1:size(final,1);
            fprintf(fid,formatspec,final{row,:});
        end

        fclose(fid);

        cd(origin)

    end
    
end

    % Done. Close screen and exit
    Screen('CloseAll');
  
end
