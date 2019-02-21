function updateContents(folder)
%UPDATECONTENTS Create a Contents.m file including subdirectories
% 
%   IOSR.GENERAL.UPDATECONTENTS scans through the current directory, and
%   its subdirectories, and builds a Contents file similar to Matlab's
%   report-generated Contents.m files. Any existing Contents.m file will be
%   overwritten.
%   
%   IOSR.GENERAL.UPDATECONTENTS(FOLDER) scans through the directory FOLDER.
% 
%   Typing
%   
%       help(FOLDER)
%   
%   or
%   
%       help path/to/folder
%   
%   will display Contents.m in the Command Window, and display links to the
%   help for any functions that are in Matlab's search path.
% 
%   NB: Do not use Matlab's Contents Report generator to edit the
%   Contents.m file. Execute this function to update it.

%   Copyright 2016 University of Surrey.

    % apply function in current directory
    if nargin<1
        folder = cd;
    end

    % check input is valid
%     assert(ischar(folder), 'iosr:updateContents:invalidPath', '''directory'' should be a charater array (string)')
%     assert(exist(folder,'dir')==7, 'iosr:updateContents:invalidPath', [folder ' does not exist'])

    % check last character in path is not filesep (e.g. '/')
    if strcmp(folder(end),filesep)
        folder = folder(1:end-1);
    end
    fIX = strfind(folder,filesep); fIX = fIX(end);

    % name of file to create
    filename = 'Contents.m';
    
    % Name of the folder
    [~,name] = fileparts(folder);

    % delete if it already exists
    if exist([folder filesep filename],'file')==2
        delete([folder filesep filename])
    end

    % get subfolders
    dirs = getContents(folder,'filter','folders','rec',true,'path','full','sort',true);
    dirs = [{folder}; dirs];

    % get files
    files = cell(0,1);
    H1_lines = cell(0,1);
    for d = 1:length(dirs)
        temp = getContents(dirs{d},'filter','files','sort',true);
        if ~isempty(temp)
            temp = temp(cellfun(@(x) isempty(strfind(x,'~')),temp)); % remove temporary files
            temp = temp(cellfun(@(x) isempty(strfind(x,'.mex')),temp)); % remove compiled mex files
            temp = temp(cellfun(@(x) ~strcmp(x,filename),temp)); % remove Contents.m
            H1_lines = [H1_lines; {''}; {''}]; %#ok<AGROW> % insert blank lines where no functions will be
            % determine package prefix
            pkgprefix = strrep(dirs{d},[filesep '+'],'.');
            pkgprefix = strrep(pkgprefix,[filesep '@'],'.');
            dots = strfind(pkgprefix,'.');
            if ~isempty(dots)
                pkgprefix = [pkgprefix(dots(1)+1:end) '.'];
            else
                pkgprefix = '';
            end
            for f = 1:length(temp) % read H1 lines
                H1_lines = [H1_lines; {get_H1_line([dirs{d} filesep temp{f}])}]; %#ok<AGROW> % add H1 lines
                % remove extension from and add package prefix to m-files
                [~,fname,ext] = fileparts(temp{f});
                if strcmpi(ext,'.m')
                    temp{f} = [pkgprefix fname];
                end
            end
            files = [files; {''}; {upper(dirs{d}(fIX+1:end))}; temp;]; %#ok<AGROW> % add filenames
        end
    end

    % longest file name (so appropriate space can be added between files and H1 lines
    longest_word = max(cellfun(@length,files(cellfun(@(x) ~isempty(x),H1_lines))));

    % write to output
    nrows = length(files);
    fid = fopen(filename, 'w'); % open file for writing
    fprintf(fid, '%s\n%% \n', ['% ' upper(name)]);
    fprintf(fid, '%s\n', ['%   Contents file for ' upper(folder(fIX+1:end)) ' and its subfolders.']);
    for row=1:nrows
        if isempty(H1_lines{row})
            fprintf(fid, '%s\n', ['%   ' files{row,:}]);
        else
            rowfilename = files{row,:};
            [~,name,ext] = fileparts(rowfilename);
            if strcmpi(ext,'.m') % remove extension from m files
                rowfilename = name;
            end
            fprintf(fid, '%s\n',['%   ' rowfilename repmat(' ',1,longest_word-length(rowfilename)) ' - ' H1_lines{row,:}]);
        end
    end
    fprintf(fid, '%%    \n%%   %s on %s at %s.\n', 'This file was generated by updateContents.m',datestr(datetime('now'),'dd mmm yyyy'),datestr(datetime('now'),'HH:MM:SS'));
    fclose(fid);

end

function H1_line = get_H1_line(filename)
%GET_H1_LINE get the H1 line for a file

    [~,name,ext] = fileparts(filename);
    H1_line = ''; % default output
    if strcmp(ext,'.m')
        fid = fopen(filename); % open file
        tline = fgetl(fid); % read first line
        while ischar(tline)
            k = strfind(tline,'%'); % find comment
            if ~isempty(k) % if it is found
                k = k(1);
                ispercents = false(size(tline(k:end)));
                ispercents(strfind(tline(k:end),'%'))=true;
                start = k+find(~(isspace(tline(k:end)) | ispercents),1,'first')-1;
                if ~isempty(start)
                    tline = tline(start:end); % remove leading space/percent
                    IX = strfind(lower(tline),lower(name));
                    if ~isempty(IX)
                        if IX(1)==1
                            tline = tline(length(name)+1:end); % remove function name
                        end
                        tline = strtrim(tline); % remove any leading/trailing space
                    end
                    H1_line = tline;
                    H1_line = strtrim(H1_line);
                    if ~isempty(H1_line)
                        if strcmp(H1_line(end),'.') % remove trailing period
                            H1_line = H1_line(1:end-1);
                        end
                        H1_line(1) = upper(H1_line(1)); % capitalize first letter
                    end
                end
                tline = -1; % set tline to numeric
            else
                tline = fgetl(fid); % read next line
            end
        end
        fclose(fid);
    end

end
