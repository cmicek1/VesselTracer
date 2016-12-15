classdef Segment
    properties
        
        % Surface data
        x
        y
        z
        
        % Shape parameter
        epsilon
        
        % Center (x, y, z)
        mu
        
        % Scale (x, y, z)
        sigma
        
        % Orientation (z, y, x)
        phi
        
        % Rotation Matrix (ZYX)
        R
        
        % Quaternion (w, x, y, z)
        q
        
        % For Jacobian calculation
        G
        
        % Background and foreground image intensities (original estimates)
        IF0
        IB0
        
        % Background and foreground image intensities (updated)
        IF
        IB
        
        % MAD statistic
        alpha
        
        % Max iterations for initialization
        initIterMax
        
        % Log-likelihood
        L
    end
    methods (Access = protected)
        % Marginal distribution for Laplacian noise
        function val = f(obj, arg)
            val = obj.alpha / 2 * exp(-obj.alpha * abs(double(arg)));
        end
        
        function J = Jacobian(obj, start, point)
            d_mu = eye(3);
            d_sigma = zeros(3, 3);
            d_q = zeros(3, 4);
            d_epsilon = zeros(3, 1);
            if ~start
                sig = calcEllipse(point(1), point(2), point(3), obj.epsilon);
                u = point' / norm(point);
                qR = quat2rotm(obj.q);
                d_sigma = qR * diag(sig * u);
                p_vec = diag(obj.sigma) * sig * u;
                P = [0, -p_vec(3), p_vec(2);
                     p_vec(3), 0, -p_vec(1);
                     -p_vec(2), p_vec(1), 0];
                d_q = -qR * P * obj.G;
                d_epsilon = qR * diag(obj.sigma) * calc_dsig(point(1),...
                    point(2), point(3), obj.epsilon) * u;
            end
            J = [d_mu, d_sigma, d_q, d_epsilon];
        end
        function obj = updateG(obj)
            obj.G = 2 * [-obj.q(2), obj.q(1), obj.q(3), -obj.q(4);
                -obj.q(3), -obj.q(4), obj.q(1), obj.q(2);
                -obj.q(4), obj.q(3), -obj.q(2), obj.q(1)];
        end
    end
    methods
        
        function obj = Segment(varargin)
            if nargin > 0
                obj.x = varargin{1};
                obj.y = varargin{2};
                obj.z = varargin{3};
            end
            obj.epsilon = 0.25;
            obj.mu = [0; 0; 0];
            obj.sigma = [3; 3; 4.5];
            obj.phi = [0; 0; 0];
            obj.R = eye(3);
            obj.q = [1, 0, 0, 0];
            obj = obj.updateG();
            obj.IF0 = 0;
            obj.IB0 = 255;
            obj.alpha = 0;
            obj.initIterMax = 20;
        end
        
        function obj = translate(obj, point, imSz)
            obj.mu = obj.mu + point;
            obj.mu(1) = obj.mu(1) * (obj.mu(1) >= 1 && obj.mu(1) <= imSz(1)) + ...
                (obj.mu(1) < 1) + (obj.mu(1) > imSz(1)) * imSz(1);
            obj.mu(2) = obj.mu(2) * (obj.mu(2) >= 1 && obj.mu(2) <= imSz(2)) + ...
                (obj.mu(2) < 1) + (obj.mu(2) > imSz(2)) * imSz(2);
            obj.mu(3) = obj.mu(3) * (obj.mu(3) >= 1 && obj.mu(3) <= imSz(3)) + ...
                (obj.mu(3) < 1) + (obj.mu(3) > imSz(3)) * imSz(3);
            obj.x = obj.x + point(1);
            obj.y = obj.y + point(2);
            obj.z = obj.z + point(3);
        end
        
        function obj = rotate(obj, angle)
            if nargin > 1
                obj.q = angle;
                obj.R = quat2rotm(angle);
                obj.phi = rotm2eul(obj.R);
            else
                obj.R = quat2rotm(obj.q);
                obj.phi = rotm2eul(obj.R);
            end
            rotated = obj.R * [obj.x; obj.y; obj.z];
            obj.x = rotated(1, :);
            obj.y = rotated(2, :);
            obj.z = rotated(3, :);
            obj = obj.updateG();
        end
        
        function obj = scale(obj, scalevec)
            if nargin > 1
                obj.sigma = scalevec;
            else
                scalevec = obj.sigma;
            end
            obj.x = obj.x * scalevec(1);
            obj.y = obj.y * scalevec(2);
            obj.z = obj.z * scalevec(3);
        end
        
        function obj = shape(obj, sz, epsilon)
            if nargin > 2
                obj.epsilon = epsilon;
            else
                epsilon = obj.epsilon;
            end
            start = ellipsoidInit(epsilon, 1);
            obj.x = start.XData;
            obj.y = start.YData;
            obj.z = start.ZData;
            clear start
            obj = obj.scale(obj.sigma).rotate(obj.q).translate(obj.mu, sz);
        end
        
        
        function obj = intensity_est(obj, image, outsideWin, start)
            minr = min(obj.x);
            minc = min(obj.y);
            minz = min(obj.z);
            
            maxr = max(obj.x);
            maxc = max(obj.y);
            maxz = max(obj.z);
            
            [numrow, numcol, numz] = size(image);
            
            
            lr = floor(max([1, minr]));
            ur = floor(min([maxr, numrow]));
            
            lc = floor(max([1, minc]));
            uc = floor(min([maxc, numcol]));
            
            lz = floor(max([1, minz]));
            uz = floor(min([maxz, numz]));
            
            % Get bounding box for ellipsoid (easier/more efficient than convex
            % hull)
            insideRegion = image(lr:ur, lc:uc, lz:uz);
            
            % Estimate foreground intensity
            if start
                obj.IF0 = median(insideRegion(:));
                obj.IF = obj.IF0;
            else
                % Use skipped median M-estimator for future estimates
                temp = insideRegion(:) - obj.IF0;
                s = 0.5 * abs(obj.IF0 - obj.IB0);
                temp = temp(temp <= s);
                obj.IF = median(temp);
            end
            
            % Calculate MAD statistic
            obj.alpha = mad(double(insideRegion(:)));
            
            lr = max([1, lr - outsideWin]);
            ur = min([ur + outsideWin, numrow]);
            
            lc = max([1, lc - outsideWin]);
            uc = min([uc + outsideWin, numcol]);
            
            lz = max([1, lz - 2]);
            uz = min([uz + 2, numz]);
            
            % Get large region outside ellipsoid
            outsideRegion = image(lr:ur, lc:uc, lz:uz);
            
            % Estimate background intensity
            if start
                obj.IB0 = median(outsideRegion(:));
                obj.IB = obj.IB0;
            else
                % Use skipped median M-estimator for future estimates
                temp = insideRegion(:) - obj.IB0;
                s = 0.5 * abs(obj.IF0 - obj.IB0);
                temp = temp(temp <= s);
                obj.IB = median(temp);
            end
        end
        
        function [delL, L] = Lgrad(obj, image, start) 
            try
                hull = convhull(obj.x, obj.y, obj.z, 'Simplify', true);
            catch EXCEPT
                delL = zeros(11, 3);
                L = NaN;
                return
            end
            set(0,'DefaultFigureVisible','off');
            t = trisurf(hull, obj.x, obj.y, obj.z);
            set(0,'DefaultFigureVisible','on');
            np = reducepatch(t, 0.2);
            [numrow, numcol, numz] = size(image);
            
            Aind = np.faces(:, 1);
            Bind = np.faces(:, 2);
            Cind = np.faces(:, 3);
            
            Av = np.vertices(Aind, :);
            Bv = np.vertices(Bind, :);
            Cv = np.vertices(Cind, :);
            
            c = (Av + Bv + Cv) / 3;
            
            % Threshold
            c(c < 1) = 1;
            tmp1 = c(:, 1);
            tmp1(tmp1 > numrow) = numrow;
            tmp2 = c(:, 2);
            tmp2(tmp2 > numcol) = numcol;
            tmp3 = c(:, 3);
            tmp3(tmp3 > numz) = numz;
            c = [tmp1, tmp2, tmp3];
            
            clear tmp1 tmp2 tmp3
            c_round = int64(floor(c));
            li = sub2ind([numrow, numcol, numz], c_round(:, 1),...
                c_round(:, 2), c_round(:, 3));
            
            IC = image(li);
            aVec = cross(Bv - Av, Cv - Av);
            lens = sqrt(sum(abs(aVec).^2,2));
            surfnorms = aVec ./ lens;
            areas = 0.5 * lens;
            
            F = log(obj.f(IC - obj.IF) ./ obj.f(IC - obj.IB));
            F(isinf(F)) = 0;
            if start
                J = permute(obj.Jacobian(true)', [1, 2, 3]);
                J = repmat(J, 1, 1, length(F));
                temp = (F .* areas .* surfnorms)';
                temp = permute(temp, [1, 3, 2]);
                delL = sum(mtimesx(J, temp, 'SPEEDOMP'), 3);
                if any(isnan(delL))
                    delL(isnan(delL)) = 0;
                end
                L = sum(F);
            else
                J = NaN(11, 3, length(F));
                for i = 1:size(c, 1)
                    J(:, :, i) = obj.Jacobian(false, c(i, :))';
                end
                temp = (F .* areas .* surfnorms)';
                temp = permute(temp, [1, 3, 2]);
                delL = sum(mtimesx(J, temp, 'SPEEDOMP'), 3);
                L = sum(F);
            end
        end
    end
end

function sigma = calcEllipse(x, y, z, epsilon)
sigma = ((x.^2 + y.^2).^(1/epsilon) + z.^(2/epsilon)).^(-epsilon/2);
end

function dsde = calc_dsig(x, y, z, eps)
dsde = 2/eps^3 * ((x.^2 + y.^2).^(1/eps) + z.^(2/eps)).^(-(eps + 2)/eps).*...
    (eps*((x.^2 + y.^2).^(1/eps) + z.^(2/eps))*...
    log((x.^2 + y.^2).^(1/eps) + z.^(2/eps)) + (x.^2 + y.^2).^(1/eps).*...
    log(x.^2 + y.^2) + 2*z.^(2/eps).*log(z));
end