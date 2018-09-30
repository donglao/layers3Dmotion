function composed = compose_flow_single( flow1,flow2)
[h, w] = size( flow1( :,:,1 ) );
[x, y] = meshgrid( 1:w, 1:h );
temp1 = interp2( x, x+flow2( :,:,1 ), y+flow2( :,:,2 ) );
temp2 = interp2( y, x+flow2( :,:,1 ), y+flow2( :,:,2 ) );
temp1 = interp2( temp1, x+flow1( :,:,1 ), y+flow1( :,:,2 ) );
temp2 = interp2( temp2, x+flow1( :,:,1 ), y+flow1( :,:,2 ) );
composed( :,:,1 )= temp1-x;
composed( :,:,2 )= temp2-y;
end