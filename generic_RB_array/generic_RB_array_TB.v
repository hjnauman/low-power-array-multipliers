`timescale 1ns / 1ps

module generic_RB_mult_TB;
  wire [15:0] p;
  reg [7:0] a;
  reg [7:0] x;
  reg [7:0] i;

  ArrayMultiplier_row_bypass_generic #(.m(8), .n(8)) am (p, a, x);
//  initial $monitor("a=%b,x=%b,p=%b,", a, x, p);

  initial
  begin
    i = 0;  
    x = 8'b11111111;
    a = 8'b00000001;
    #10 x = 8'b11111110;
    #10 x = 8'b11111101;
    #10 x = 8'b11111011;
    #10 x = 8'b11110111;
    #10 x = 8'b11101111;
    #10 x = 8'b11011111;
    #10 x = 8'b10111111;
    #10 x = 8'b01111111;
    #10 $finish;

//    for ( i=0 ; i< 2'd16; i = i + 1) begin
//        #10 a = a + 1'b1;
//    end
//    #1 $finish;
    end
endmodule