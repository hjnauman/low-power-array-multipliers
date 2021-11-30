`timescale 1ns / 1ps //check to see that it's coming through
module Full_Adder(input a, b, carry_in, output reg sum, carry_out);
    always @(a or b or carry_in)
    begin
        sum = a ^ b ^ carry_in;
        carry_out = a & b | (a^b) & carry_in;
    end
endmodule

module mux(input a, b, select, output reg out);
    always @(a or b or select)
    begin
        out = (a&~select)|(b&select);
    end
endmodule 

module tri_state(input a, enable, output reg out);
    always @(a or enable)
    begin
    if (enable == 0)
        out = 1'bz;        
    else
        out = a;
    end
endmodule

module two_in_and_gate(input a, b, output reg z);
    always @(a or b)
    begin
        z = a & b;
    end
endmodule

module FACB_four_input(input a_old, a_new, b_old, b_new, enable, carry_in, output carry_out, sum);
    wire FA_in0, FA_in1, FA_out, tri_in0, tri_in1; //, Mux0_in, Mux1_in;
    two_in_and_gate and0(.a(a_old), .b(b_new), .z(tri_in0));
    two_in_and_gate and1(.a(a_new), .b(b_old), .z(tri_in1));
    tri_state tri_0(.a(tri_in0), .enable(enable), .out(FA_in0));
    tri_state tri_1(.a(tri_in1), .enable(enable), .out(FA_in1));
    Full_Adder FA(.a(FA_in0), .b(FA_in1), .carry_in(carry_in), .sum(FA_out), .carry_out(carry_out));
    mux mux0(.a(tri_in1), .b(FA_out), .select(enable), .out(sum));
endmodule

module FACB_three_input(input a_old, b_new, in1, enable, carry_in, output carry_out, sum);
    wire FA_in0, FA_in1, FA_out, tri_in0; //, Mux0_in, Mux1_in;
    two_in_and_gate and0(.a(a_old), .b(b_new), .z(tri_in0));
    tri_state tri_0(.a(tri_in0), .enable(enable), .out(FA_in0));
    tri_state tri_1(.a(in1), .enable(enable), .out(FA_in1));
    Full_Adder FA(.a(FA_in0), .b(FA_in1), .carry_in(carry_in), .sum(FA_out), .carry_out(carry_out));
    mux mux0(.a(in1), .b(FA_out), .select(enable), .out(sum));
endmodule

module Full_Adder_LR_three_in(input a, and_value, b, carry_in, output reg sum, carry_out);
    always @(a or b or and_value or carry_in)
    begin
        sum = a ^ (b&and_value) ^ carry_in;
        carry_out = a & (b&and_value) | (a^(b&and_value)) & carry_in;
    end
endmodule

module Full_Adder_LR_four_in(input in0, in1, in2, in3, carry_in, output reg sum, carry_out);
    always @(in0 or in1 or in2 or in3 or carry_in)
    begin
        sum = (in0&in1) ^ (in2&in3) ^ carry_in;
        carry_out = (in0&in1) & (in2&in3) | ((in0&in1)^(in2&in3)) & carry_in;
    end
endmodule

module ArrayMultiplier_column_bypass_generic(product, a, b);
    parameter m = 4;
    parameter n = 4;
    output wire [m+n-1:0] product;
    input[m-1:0] a;
    input[n-1:0] b;

    wire c_partial[(m-1)*(n)-1:0];
    wire s_partial[(m-1)*(n)-1:0];

    //generate first product bit
    two_in_and_gate two_and(.a(a[0]), .b(b[0]), .z(product[0]));

    //generate first row of the multiplier
    genvar i;
    generate
        for(i=0; i<m-1; i=i+1)
        begin
            //FACB_four_input(input a_old, a_new, b_old, b_new, enable, carry_in, output carry_out, sum);
            FACB_four_input FA_first(.a_old(a[i]), .a_new(a[i+1]), .b_old(b[0]), .b_new(b[1]), .enable(a[i]), .carry_in(1'b0),
             .sum(s_partial[i]), .carry_out(c_partial[i]));
        end
    endgenerate

    //generate middle rows of the multiplier except last column
    genvar j, k;
    generate
        for(k=0; k<n-2; k=k+1) //row
        begin
            for(j=0; j<m-2; j=j+1)  //column
            begin
                //FACB_three_input(input a_old, b_new, in1, enable, carry_in, output carry_out, sum);
                FACB_three_input FA_middle(.a_old(a[j]), .b_new(b[k+2]), .in1(s_partial[((m-1)*k)+(j+1)]), .enable(a[j]), .carry_in(c_partial[((m-1)*k)+j]),
                 .sum(s_partial[((m-1)*(k+1))+j]), .carry_out(c_partial[((m-1)*(k+1))+j]));
            end
        end
    endgenerate

    //generate middle lines of the multiplier - last column of each row
    genvar z;
    generate
        for(z=0; z<n-2; z=z+1)
        begin
            //FACB_four_input(input a_old, a_new, b_old, b_new, enable, carry_in, output carry_out, sum);
            FACB_four_input FA_middle_last_column(.a_old(a[m-2]), .a_new(a[m-1]), .b_old(b[z+1]), .b_new(b[z+2]), .enable(a[m-2]), .carry_in(c_partial[(m-1)*(z)+(m-2)]),
             .sum(s_partial[((m-1)*(z+1))+(m-2)]), .carry_out(c_partial[((m-1)*(z+1))+(m-2)]));
        end
    endgenerate

    //generate last line of the multiplier
    genvar l;
    generate
        for(l=0; l<n-1; l=l+1)
        begin
            if(l==0) //first column
            begin
                //Full_Adder_LR_three_in(input a, and_value, b, carry_in, output reg sum, carry_out);
                Full_Adder_LR_three_in Last_row_first_column(.a(s_partial[((m-1)*(n-2))+1]), .and_value(a[l]), .b(c_partial[((m-1)*(n-2))+l]), .carry_in(1'b0),
                 .sum(s_partial[((m-1)*(n-1))+l]), .carry_out(c_partial[((m-1)*(n-1))+l]));
            end else if (l < m-2)
            begin
                Full_Adder_LR_three_in Last_row_middle_columns(.a(s_partial[((m-1)*(n-2))+l+1]), .and_value(a[l]), .b(c_partial[((m-1)*(n-2))+l]), .carry_in(c_partial[((m-1)*(n-1))+l-1]),
                 .sum(s_partial[((m-1)*(n-1))+l]), .carry_out(c_partial[((m-1)*(n-1))+l]));
            end else begin
                //Full_Adder_LR_four_in(input in0, in1, in2, in3, carry_in, output reg sum, carry_out);
                Full_Adder_LR_four_in Last_row_last_column(.in0(a[m-1]), .in1(b[n-1]), .in2(a[m-2]), .in3(c_partial[((m-1)*(n-2))+l]), .carry_in(c_partial[((m-1)*(n-1))+l-1]), 
                 .sum(s_partial[((m-1)*(n-1))+l]), .carry_out(c_partial[((m-1)*(n-1))+l]));
            end
        end
    endgenerate

    //end product bits from first and middle cells
    generate
        for(i=0; i<n-1; i=i+1)
        begin
            buf (product[i+1], s_partial[(m-1)*i]);
        end
    endgenerate

    //end product bits for last line of cells
    generate
        for(i=0; i<n-1; i=i+1)
        begin
            buf (product[i+n], s_partial[((m-1)*(n-1))+i]);
        end
    endgenerate

    //msb of product
    buf (product[m+n-1], c_partial[((m-1)*(n))-1]);
endmodule