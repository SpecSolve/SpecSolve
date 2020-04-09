function [percent, elapsed] = parfor_progress(N, varargin)
%PARFOR_PROGRESS Progress and duration for parallel executions (e.g. parfor)
%
% PARFOR_PROGRESS works by storing progress and start time in a temporary
% file. This is needed because parallel threads cannot communicate with one
% another so there is no simple way to know the execution order. This
% function is safe to use for large number of iterations as only 3 numbers
% are written to the file, i.e. it stays small (read and written quickly).
% Obviously, it works for non-parallel iterations too.
%
% FILE = PARFOR_PROGRESS(N) initializes the progress monitor for a set
%   of N upcoming iterations, parallel or not.
%
% PARFOR_PROGRESS(FILE) updates the progress inside your parfor loop and
%   displays an updated progress bar and estimated remaining time.
%
% DELETE(FILE) deletes the parfor progress file. A better way is to set a 
%   cleanup code C=ONCLEANUP(@() DELETE(FILE)) instead, in case an error or
%   user-break occurs during parfor execution, or else files may stack up.
%
% To suppress output just ask for a return value, e.g.
%   PERCENT = PARFOR_PROGRESS(FILE) or
%   [PERCENT, ELAPSED] = PARFOR_PROGRESS(FILE) which returns the progress
%   fraction and elapsed time in seconds.
%
% EXAMPLE:
%
%   N = 100;
%   pf = parfor_progress(N);
%   fin = onCleanup(@() delete(pf));    % delete file on exit/error/break
%   parfor i=1:N
%       pause(rand);                    % Replace with real code
%       parfor_progress(pf);
%   end
%
% See also PARFOR.
%
% Author: Bogdan Roman, http://www.damtp.cam.ac.uk/user/abr28, 2014
% inspired by the work of Jeremy Scheff, http://www.jeremyscheff.com

narginchk(1, 2);

if isnumeric(N) && N > 0
    if nargin > 1
        file = varargin{1};
    else
        file = [tempname '_parfor.txt'];
    end
    fid = fopen(file, 'w');
    if fid < 0, error('Could not open file for writing (perms?): %s', file); end
    % write N, start time (0.1s resolution), and iteration (0 for now)
    progress = [N floor(now*864000) 0];
    fprintf(fid, '%d\n%d\n0\n', progress);
    fclose(fid);
    computeprogress(progress, false, true);
    percent = file;
    elapsed = 0;
else
    file = N;
    if ~ischar(file) || ~exist(file, 'file')
        error('Not initialized. See HELP PARFOR_PROGRESS.');
    end
    % update (read and write) in one go
    fid = fopen(file, 'r+');
    progress = fscanf(fid, '%f');                   % read the 3 values
    progress(3) = progress(3) + 1;                  % update iteration number
    fseek(fid, 0, 'bof');
    fprintf(fid, '%d\n%d\n%d\n', round(progress));  % write back to file
    fclose(fid);
    [percent, elapsed] = computeprogress(progress, true, nargout == 0);
end

end

%---------------- PRIVATE FUNCTIONS ----------------%

function [percent, elapsed] = computeprogress (progress, update, show)
    elapsed = (now - progress(2)/864000) * 86400;   % compute elapsed seconds
    percent = progress(3) / progress(1);
    if percent == 0
        duration = 0;
        remaining = 0;
    else
        duration = elapsed / percent; % TODO: improve this crude estimate, exp smoothing filter?
        remaining = duration - elapsed;
    end
    if show
        r = humantime(remaining);
        e = humantime(elapsed);
        t = humantime(duration);
        s = sprintf('%8.2f%%, %s (el), %s (rem), %s (tot)\n', ...
                percent * 100, e, r, t);
        if update, back = repmat(char(8),1,length(s)); else back = ''; end
        fprintf('%s%s', back, s);
%         % if not GUI mode then show one line per tick (useful for log files)
%         %screenSize = get(0,'ScreenSize');
%         %if isequal(screenSize(3:4),[1 1])
%         if ~usejava('jvm') || ~feature('ShowFigureWindows')
%             fprintf('%6.2f%%, %s (el), %s (rem), %s (tot)\n', ...
%                 percent * 100, e, r, t);
%         else
%             width = 50; % width of progress bar
%             ticks = round(percent*width);
%             if update, back = repmat(char(8), 1, (width+8+length(r)+length(e)+4)); else back = ''; end
%             disp([back, sprintf('%3d%%',round(percent*100)), ' [', ...
%                repmat('=', 1, ticks), repmat(' ', 1, width - ticks), '] ' e ' / ' r]);
%         end
    end
end

function t = humantime (s)
    if s < 60, t = sprintf('%4.1fs', s);
    elseif s < 3600, t = sprintf('%4.1fm', s/60);
    elseif s < 86400, t = sprintf('%4.1fh', s/3600);
    elseif s < 604800, t = sprintf('%4.1fd', s/86400);      % 86400 = 1 day = 24 * 3600
    elseif s < 2629800, t = sprintf('%4.1fw', s/604800);    % 604800 = 1 week = 7 * 86400
    elseif s < 31557600, t = sprintf('%4.1fM', s/2629800);  % 2629800 = 1 average month = 365.25/12 * 86400
    else t = sprintf('%4.1fy', s/31557600);                 % 31557600 = 1 average year = 365.25 * 86400
    end
    %t = sprintf('%3dh%02dm%02ds', floor(s/3600), mod(floor(s/60),60), mod(floor(s),60));
end
