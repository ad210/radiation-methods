arr = csvread('my_spectrum_all_000.csv');

imagesc(arr(:,2:1024))
colorbar

max(max(arr(:,2:1024))) * 1E30