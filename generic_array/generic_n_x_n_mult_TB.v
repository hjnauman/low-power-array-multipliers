`timescale 1ns / 1ps

module generic_n_x_n_mult_TB;
  wire [15:0] p;
  reg [7:0] a;
  reg [7:0] x;
  reg [7:0] i, j;

  ArrayMultiplier_generic #(.m(8), .n(8)) am (p, a, x);
//  initial $monitor("a=%b,x=%b,p=%b,", a, x, p);

  initial
  begin
    a = 0;
    x = 0;
    i = 0;
    j = 0;
    for ( i=0 ; i< 64; i = i + 1) begin
          a = a + 1'b1;
        x = 8'b00000000;
        for( j = 0; j<64; j = j + 1) begin
            #2 x = x + 1'b1;
        end
    end
    #10 $finish;
    end
endmodule