function [varargout] = CalcGrainSize(varargin)
    
    % Create Grain (and Twin--if necessary) structure for data output
    E2627 = [];
    Twins = [];
    
    % First object is always the ebsd
    ebsd = varargin{1};
    % Second object is always the phase to extract
    PhaseId = varargin{2};
    % Third object is always the flag of whether to merge the twins into 
    % their parent orientation
    MergeTwins = varargin{3};
    
    % Fourth object would be the grain boundary misorientation tolerance
    % angle (what classifies the GB interface: default is 5 degrees)
    if length(varargin) == 4
        GB_Tol = varargin{4} * degree;
    else
        GB_Tol = 5*degree;
    end
    
    PxTol = 100;
    % NOTE 1: Pixel tolerance can also be input argument, although default
    % ASTM standards is 100 pixels;
    
    % Determine number of phases present in ebsd dataset. 
    EbPh = ebsd(ismember(ebsd.phaseId,PhaseId));
    
    % Compute grain boundaries and characteristics for Phase and
    % extract all grain ids for Phase
    if MergeTwins
        [GrnPh,EbPh,S3,S9] = find_twin_gbs(EbPh, PhaseId);
        Twins.S3 = S3;
        Twins.S9 = S9;
    else
        [GrnPh,EbPh.grainId] = calcGrains(EbPh,'unitCell','angle',GB_Tol);
    end
    GrnId = GrnPh.id;

    % Initiate indexing for which grains fail
    GrnFail = (GrnPh.grainSize < PxTol);

    % Ebsd step size (ASTM scaling factors)
    stepsz = EbPh.scanUnit;
    if stepsz == 'um'
        delta = 1e-6;
        scale = 1e6;
    elseif stepsz == 'mm'
        delta = 1e-3;
        scale = 1;
    else
        warning('Ebsd step size not recorded')
    end

    % Extract the unique grain IDs bordering the IPF map and ignore
    % from calculations (along with smaller grains)
    GrnBB = GrnPh.isBoundary;
    GrnFail(GrnBB) = 1;

    % Also extract grains included in ASTM
    GrnPass = GrnId(logical(abs(GrnFail-1)));

    % Extract area and % multiply measured grain area by (step size)^2
    ASTM_GrnArea = area(GrnPh(GrnPass)) * delta^2;

    % If hex-grid, modify to follow ASTM E2627 standards
     if length(EbPh.unitCell) == 6
         ASTM_GrnArea = ASTM_GrnArea * sqrt(3)/2;
     end

    % Compute ASTM E2627 grain size number
    ASTM_E2627_G = -3.3223*log10(mean(ASTM_GrnArea*scale))-2.995;

    % Collect statistics from grain size calculation
    Sz = size(GrnPh(GrnPass),1);
    freqGS = 1/Sz * sum(ASTM_GrnArea);
    StDev = sqrt(1/(Sz-1) * sum((ASTM_GrnArea - freqGS).^2));
    StError = StDev / sqrt(Sz);
    CI95 = tstat3(max(Sz-1,1), 0.975, 'inv');
    yCi95 = bsxfun(@times, StError, CI95(:));
    
    % Store ASTM E2627 Grain Size Calculations
    E2627.AcceptableGrainIds = GrnPass;
    E2627.GrainSize          = ASTM_E2627_G;
    E2627.Ci95               = yCi95;
    
    % Calculate the ASTM E112_13 grain size number
    [E112_13] = Calculate_ASTM_E112_13(GrnPh,EbPh);
    
    % Pack the kids up and go home
    varargout{1} = EbPh;
    varargout{2} = GrnPh;
    varargout{3} = E2627;
    varargout{4}  = E112_13;
    if isempty(Twins) == 0
        varargout{5} = Twins;
    end

end

function [E112_13] = Calculate_ASTM_E112_13(varargin)

%%                      Function Description

% This function follows ASTM E112-13 standards (for grains) by combining 
% annealling twins with the parent s.t. the parent-twin system becomes a 
% single entity whose grain boundaries are defined by the parent grain.

% Single intercept lines are performed on any packet, block, and sub-block 
% boundaries for further analysis if the material is steel by randomly
% assigning a start and stop on the x_min and y_max edge, and then
% performing the same sort of standard segmentation analysis.

% The grain boundary lineal intercept method involves 12 total lineal 
% intercepts, with 8 vertical lines spanning the x-axis coordinates and 
% 4 horizontal lines spanning the y-axis, all lines being equally 
% distributed from each other.


%%                      Function Implementation

    % Grab grains and ebsd
    grains = varargin{1};
    ebsd   = varargin{2};

    % Assign myEBSD and Grain structures (called Feature here)
    Bounds = grains.boundary;
    
    % If user wants to define their own line segment, use their
    % coordinates. Else, randomly confine a line to the minimum x value 
    % (left edge), a y value contained within this edge, and two randomly 
    % chosen x,y endpoints.
    if length(varargin)==3
        X1 = varargin{3}(1);
        X2 = varargin{3}(2);
        Y1 = varargin{3}(3);
        Y2 = varargin{3}(4);

        % Establish the line endpoints
        L1 = [X1,Y1];
        L2 = [X2,Y2];
        % Length of line
        lenLine = sqrt((L2(1)-L1(1))^2+(L2(2)-L1(2))^2);
        
        % Determine number of intersections and line-segment lengths
        [X,Y,SortedLineSegs] = Bounds.intersect(L1,L2);
        Xint = X(find(~isnan(X)));
        Yint = Y(find(~isnan(Y)));
        % Number of intersections
        P_i = length(SortedLineSegs);
        
        % Average grain size
        li = lenLine / P_i;
        % Surface Area
        Sv = 2/lenLine*1000;
        
        % Calculate ASTM Grain Size Number Per ASTM E112-13
        ASTM_Num = -6.643856*log10(lenLine/1000)-3.288;
        
    else % Do multiple lineal intercepts for ASTM grain size measurement
        
        % Extract triple point coordinates and set pixelated tolerance
        tpX = Bounds.triplePoints.x;
        tpY = Bounds.triplePoints.y;
        tpTol = 1;

        % Range of x and y values and lineal intercept line length
        Lx = range(ebsd.x);
        Ly = range(ebsd.y);
        TotLine = (Lx*4)+(Ly*8);

        % X-based values
        XPart = Lx / 10;
        XStrt = min(ebsd.x);
        XEnd = max(ebsd.x);   

        % Y-based values
        YPart = Ly / 6;
        YStrt = min(ebsd.y);
        YEnd = max(ebsd.y);

        % X-dependent line that runs across the x-axis coordinates at 
        % constant y values
        Xdep_X1 = ones(1,4)*XStrt-10;
        Xdep_X2 = ones(1,4)*XEnd+10;
        Xdep_Y1 = linspace(YStrt+YPart,YEnd-YPart,4);
        Xdep_Y2 = Xdep_Y1;

        % Y-dependent line that runs across the y-axis coordinates at 
        % constant x values
        Ydep_X1 = linspace(XStrt+XPart,XEnd-XPart,8);
        Ydep_X2 = Ydep_X1;
        Ydep_Y1 = ones(1,8)*YStrt-10;
        Ydep_Y2 = ones(1,8)*YEnd+10;

        % Number of lineal intercepts
        Sz = 12;

        % Preallocate arrays
        li        = zeros(Sz,1); % Line segment length
        NSegs     = zeros(Sz,1); % Number of intersections
        Sv        = zeros(Sz,1); % Surface Area
        ASTM_Nums = zeros(Sz,1); % ASTM E112-13 grain size number
        
        % Loop through all line segments and calculate grain size based on
        % lineal intercept method
        for ii = 1:Sz

            if ii > 4
                X1 = Ydep_X1(ii-4);
                X2 = Ydep_X2(ii-4);
                Y1 = Ydep_Y1(ii-4);
                Y2 = Ydep_Y2(ii-4);
            else
                X1 = Xdep_X1(ii);
                X2 = Xdep_X2(ii);
                Y1 = Xdep_Y1(ii);
                Y2 = Xdep_Y2(ii);
            end

            % Establish the line endpoints
            L1 = [X1,Y1];
            L2 = [X2,Y2];

            % Determine number of intersections and line-segment lengths
            [X,Y,~] = Bounds.intersect(L2,L1);

            % Extract the intersections
            Ints = ~isnan(X);
            % Id grains of the intersections. 
            IntGrns = Bounds.grainId(Ints,:);
            % NOTE 2: The first and last entries will always be 0 and an
            % existing grain id. These first 0's are the edges of the 
            % boundaries and count for 0.5 intersections

            % Determine first whether any repeated grain intersection Ids
            % coincide with the ends of the micrograph
            RepeatGrnEnds = zeros(length(IntGrns),1,'logical');
            EndIds = [1,length(IntGrns)];
            % Find whether any "middle" grain intersection pattern repeats
            % one of the ends of the micrograph ids.
            MiddleGrnRepeat = ismember(IntGrns(2:EndIds(2)-1,:),IntGrns(EndIds,:),'row');
            RepeatGrnEnds(2:EndIds(2)-1) = MiddleGrnRepeat;
            % NOTE 3: If this does occur, it means that either an unindexed
            % region or different phase grain is sandwiched between one of
            % end grains and another grain. We need to delete these grains 
            % and the corresponding line segments for accurate line lengths

            % Delete those grains that repeat at the end
            IntGrns(RepeatGrnEnds,:) = [];

            % Find X and Y coordinates of intersections
            Xint = X(Ints);
            Yint = Y(Ints);

            % Extract coordinates of the unique grains at intersections
            [UnqIntGrns,UnqIntGrnId,UnqIntRepIds] = unique(IntGrns,'rows','stable');
            UnqXint = Xint(UnqIntGrnId);
            UnqYint = Yint(UnqIntGrnId);

            % Look for tangents (where line runs across 2 grain boundaries
            % and intersects multiple times)
            RepCnts = histcounts(UnqIntRepIds);
            RepGrns = UnqIntGrns(RepCnts > 2,:);
            RepGrnsRev = fliplr(RepGrns);
            % If the reverse case also exists [i.e. (i,j) and (j,i), 
            % # delete (j,i) from the unique grain intersection list
            if any(ismember(RepGrnsRev,UnqIntGrns,'rows'))
                DelGrnRow = ismember(UnqIntGrns,RepGrnsRev,'rows');
                UnqIntGrns(DelGrnRow,:) = [];
            end

            % Determine the line segment lengths. Pick from the second
            % column because it ignores all cases where a grain intersects
            % with an unindexed or different phase grain
            UnqGrnSegs = unique(UnqIntGrns(:,2),'stable');
            lenline = 0;
            for  IntCnt = 1:length(UnqGrnSegs)
                [~,~,lenseg] = grains(UnqGrnSegs(IntCnt)).boundary.intersect(L2,L1);
                if length(lenseg) > 1
                    lenseg = sum(lenseg);
                end
                lenline = lenline+lenseg;
            end

            % Those flagged as tangents will now reduce the intersection
            % count by 0.5x tangent grain boundaries
            TanFlgs = ismember(UnqGrnSegs,RepGrns);
            % NOTE 4: ASTM counts tangent boundaries as 1 intersection,
            % whereas MTEX will flag as 2, so we flag here and then
            % subtract 0.5 intersections for eacht time this occurs. Ergo,
            % grain i and j show a tangent boundary, initially this would
            % count for 2 intersections, but each grain will get flagged
            % and thus 2-0.5*2 = 1 intersection.
            % NOTE 5: ASTM doesn't specify, but here, if a grain is tangent
            % with an non-phase grain boundary, that intersection only
            % counts as 0.5 total

            % Sort based on X intercepts unless x is constant, then sort by Y
            if X2 - X1 == 0
                [~,SortIntGrnId] = sort(UnqYint,'descend');
            else
                [~,SortIntGrnId] = sort(UnqXint);
            end

            % Extract the sorted, unique grain ids of intercepts and x and y 
            % coordinate intercepts. Although sorting is unneccessary, it 
            % allows us to follow the line from grain to grain either
            % horizontally or vertically
            SortedXint = UnqXint(SortIntGrnId)';
            SortedYint = UnqYint(SortIntGrnId)';

            % Number of intersections. A full intersection consists of a
            % line segment passing through two same-phase grains; if a 0
            % exists, it's either a boundary grain or a different phase
            % (unindexed). These intersections count as 0.5 intersections.
            P_i = length(UnqIntGrns) / 2 - 0.5*sum(TanFlgs);

            % Check for triple points
            tpDsts = abs((L2(2)-L1(2)).*tpX-(L2(1)-L1(2)).*tpY+L2(1)*L1(2)-L2(2)*L1(1))/...
                sqrt((L2(2)-L1(2))^2+(L2(1)-L1(1))^2);
            
            % For square pixels, it's simple
            if length(ebsd.unitCell) == 4
                % Calculate the tolerance for triple point intersection
                % distances
                tolDst = 2*tpTol*abs(ebsd.unitCell(1,1));
            else
            % For hexagonal pixels, keep it simple and only use the
            % x-distance between pixel centroids
            XminHex = min(ebsd.unitCell(:,1));
            XmaxHex = max(ebsd.unitCell(:,1));
            tolDst = tpTol*(XmaxHex-XminHex);
            end
            
            % Find the length of the number of triple points whose
            % tolerance is intersected by the lineal intercept line
            tpSegs = length(find(tpDsts < tolDst));
            
            % Adjust the number of intersections based on the number of
            % triple point intersections
            P_i = P_i + tpSegs*0.5;
              
            % Now check for any lineal intercepts that directly intersect a
            % triple point.
            TPtChk = zeros(length(SortedXint),1);
            % NOTE 6: The way MTEX searches for grain boundaries, triple
            % points don't necessarily (nor commonly) land on the vertices
            % of indexed pixels, but rather at points within the pixels.
            % NOTE 7: Because of this, if a triple point does happen to
            % land on the vertex of a pixel, it results in 3 detected 
            % grain boundary segmentations. Since ASTM standards call for 
            % 1.5 segmentations at triple points, we need to adjust P_i by 
            % subtracting 2 (since we already added the 0.5 in the previous
            % line)
            
            % Loop through all of the intersected coordinates and find
            % those that happen to align with a triple point coordinate
            for kk = 1:length(SortedXint)
                TPtChk(kk) = any(tpX == SortedXint(kk) & tpY == SortedYint(kk));

            end

          % If any triple points exist, modify our intersection count
          % by reducing by 2xn, where n is the number of triple points
          % intersected
            if length(find(TPtChk)) > 0
                P_i = P_i - length(find(TPtChk == 1));
            end
            
            % Mu line segment length
            if length(P_i) > 1
                keyboard
            end
            li(ii) = lenline / P_i;
            % Number of intersections
            NSegs(ii) = P_i;
            % Surface Area
            Sv(ii) = 2/li(ii)*1000;
            
            % ASTM E112-13 grain size number
            ASTM_Nums(ii) = -6.643856*log10(li(ii)/1000)-3.288;
        end
        
        % CI Functions
        CI95 = tstat3(max(Sz-1,1), 0.975, 'inv');

        % Lineal Intercept Data
        E112_13.LinealIntercept.Mu = TotLine / sum(NSegs);
        E112_13.LinealIntercept.LineLength = TotLine;
        E112_13.LinealIntercept.Sigma = std(li);
        E112_13.LinealIntercept.StErr = std(li) / sqrt(Sz);
        E112_13.LinealIntercept.CI95 =...
            bsxfun(@times, std(li)/sqrt(Sz), CI95(:));
        
        % ASTM Grain Size Data
        E112_13.GrainSize.Nums = ASTM_Nums;
        E112_13.GrainSize.Mu = mean(ASTM_Nums);
        E112_13.GrainSize.StDev = std(ASTM_Nums);
        E112_13.GrainSize.ASTM_StError = std(ASTM_Nums) / sqrt(Sz);
        E112_13.GrainSize.ASTM_CI95 =...
            bsxfun(@times, std(ASTM_Nums)/sqrt(Sz), CI95(:));
        
        % Surface Volume Data
        E112_13.SurfaceVolume.Mu = mean(Sv);
        E112_13.SurfaceVolume.Sigma = std(Sv);
        E112_13.SurfaceVolume.StError = std(ASTM_Nums) / sqrt(Sz);
        E112_13.SurfaceVolume.CI95 =...
            bsxfun(@times, std(ASTM_Nums)/sqrt(Sz), CI95(:));
        
        % Save the number of segmentations
        E112_13.NumberOfSegmentations = NSegs;
        
    end
end        

function pt = tstat3( v, tp, stat )
    % TSTAT3 computes one of three t_statistics: one-sided t-probability, or
    %        two-sided t-probability, or the inverse t-statistic in any single 
    %        call.  It does not take vectors as arguments.  
    % 
    %   INPUT ARGUMENTS: 
    %       v       Degrees-of freedom (integer)
    % 
    %       tp      T-statistic: to produce t-probability
    %                   OR
    %               Probability (alpha): to produce t-statistic
    % 
    %       stat    Desired test ? 
    %                   'one'   One-tailed t-probability
    %                   'two'   Two-tailed t-probability
    %                   'inv'   Inverse t-test
    % 
    %   OUTPUT ARGUMENT:
    %       pt      T-probability OR T-statistic, as requested in ?stat?
    % 
    %   USE:
    %       FIND ONE-TAILED PROBABILITY GIVEN t-STATISTIC & DEGREES-OF-FREEDOM
    %           p = tstat3(v, t, 'one')
    % 
    %       FIND TWO-TAILED PROBABILITY GIVEN t-STATISTIC & DEGREES-OF-FREEDOM
    %           p = tstat3(v, t, 'two')
    % 
    %       FIND ONE-TAILED t-STATISTIC GIVEN PROBABILITY & DEGREES-OF-FREEDOM
    %           t = tstat3(v, p, 'inv')
    %       
    %  
    %  
    % Star Strider ? 2016 01 24 ? 
    % % % % T-DISTRIBUTIONS ? 
    % Variables: 
    % t: t-statistic
    % v: degrees of freedom
    tdist2T = @(t,v) (1-betainc(v/(v+t^2), v/2, 0.5));                              % 2-tailed t-distribution
    tdist1T = @(t,v) 1-(1-tdist2T(t,v))/2;                                          % 1-tailed t-distribution
    % This calculates the inverse t-distribution (parameters given the
    %   probability ?alpha? and degrees of freedom ?v?: 
    t_inv = @(alpha,v) fzero(@(tval) (max(alpha,(1-alpha)) - tdist1T(tval,v)), 5);  % T-Statistic Given Probability ?alpha? & Degrees-Of-Freedom ?v?
    statcell = {'one' 'two' 'inv'};                                                 % Available Options
    nc = cellfun(@(x)~isempty(x), regexp(statcell, stat));                          % Logical Match Array
    n = find(nc);                                                                   % Convert ?nc? To Integer
    if (length(v) > 1) || (length(tp) > 1)
        error('                    ?> TSTAT3 does not take vectorised inputs.')
    elseif isempty(n)                                                                   % Error Check ?if? Block
        error('                    ?> The third argument must be either ''one'', ''two'', or ''inv''.')
    elseif (n == 3) && ((tp < 0) || (tp > 1))
        error('                    ?> The probability for ''inv'' must be between 0 and 1.')
    elseif (isempty(v) || (v <= 0))
        error('                    ?> The degrees-of-freedom (''v'') must be > 0.')
    end
    switch n                                                                        % Calculate Requested Statistics
        case 1
            pt = tdist1T(tp, v);
        case 2
            pt = tdist2T(tp, v);
        case 3
            pt = t_inv(tp, v);
        otherwise
            pt = NaN;
    end
end
