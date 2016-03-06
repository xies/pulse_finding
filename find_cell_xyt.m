function c = find_cell_xyt(pulse,x,y,frame)

assert(numel(pulse) == 1);

cells = pulse.cells;
cx = cat(2,cells.centroid_x); cy = cat(2,cells.centroid_y);
d = sqrt( (cx(frame,:) - x).^2 + (cy(frame,:) - y).^2 );
[~,I] = sort(d);

found = 0; i = 1;
while found
    found = inpolygon(x,y,c(I(i)).vertex_x,c(I(i)).vertex_y);
end
c = cells(I(i));

end