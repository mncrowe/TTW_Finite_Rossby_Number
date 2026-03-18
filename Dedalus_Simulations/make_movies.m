% creates movies from output data files, requires ffmpeg to be installed
% and 'save_frames.m' to be added to the matlab path

cd data

files = dir('snapshots*');

for file = files'

    disp(file.name)

    dataname = [file.name '/' file.name '_s1.h5'];

    datafile = h5info(dataname);
    x_name = datafile.Groups(1).Datasets(7).Name;
    x = h5read(dataname, ['/scales/' x_name]);
    Lx = -2*x(1);

    b = permute(h5read(dataname,'/tasks/b'),[2 1 3]) + 2/Lx*x;

    if ~exist([file.name '/frames'], "dir")     % create frames directory if missing
        mkdir([file.name '/frames'])
    end

    system(['rm ' file.name '/*.mp4 >/dev/null 2>&1']);  % delete old movie to avoid y/n from ffmpeg

    save_frames(b, [file.name '/frames/frame'], 'png', cmap2(), [-1 1])

    S = size(b);
    N = num2str(length(num2str(S(3))));

    system(['ffmpeg -framerate 10 -i ' file.name '/frames/frame_%0' N 'd.png ' file.name '/mov.mp4 >/dev/null 2>&1']);

    system(['rm ' file.name '/frames/*']);
    rmdir([file.name '/frames'])

end