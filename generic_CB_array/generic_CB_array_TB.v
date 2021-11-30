`timescale 1ns / 1ps

module generic_n_x_n_mult_TB;
  wire [7:0] p;
  reg [3:0] a;
  reg [3:0] x;
  reg [3:0] i, j;

  ArrayMultiplier_column_bypass_generic #(.m(4), .n(4)) am (p, a, x);
//  initial $monitor("a=%b,x=%b,p=%b,", a, x, p);

  initial
  begin
    a = 4'b1111;
    x = 4'b1111;
    #10 a = 4'b1110;
    x = 4'b1111;
    #10 a = 4'b1101;
    x = 4'b1111;
    #10 a = 4'b1011;
    x = 4'b1111;
    #10 a = 4'b0111;
    x = 4'b1111;
    #10 $finish;
//    a = 0;
//    x = 0;
//    i = 0;
//    j = 0;
//    for ( i=0 ; i< 64; i = i + 1) begin
//          a = a + 1'b1;
//        x = 8'b00000000;
//        for( j = 0; j<64; j = j + 1) begin
//            #2 x = x + 1'b1;
//        end
//    end
//    #10 $finish;
    end
endmodule