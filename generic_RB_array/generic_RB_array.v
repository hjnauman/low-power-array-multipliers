`timescale 1ns / 1ps

module Full_Adder(input a, b, carry_in, output reg sum, carry_out);
    always @(a or b or carry_in)
    begin
        sum = a ^ b ^ carry_in;
        carry_out = a & b | (a^b) & carry_in;
    end
endmodule

module Full_Adder_Last(input a_row, b_row, b, carry_in, output reg sum, carry_out);
    always @(a_row or b_row or b or carry_in)
    begin
        sum = (a_row&b_row) ^ b ^ carry_in;
        carry_out = (a_row&b_row) & b | ((a_row&b_row)^b) & carry_in;
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

module Full_Adder_RB_Four_In(input last_carry, a_old, a_new, b_old, b_new, carry_in, output carry_out, sum);
    wire FA_in0, FA_in1, FA_carry_in, FA_out0, FA_out1, and0_out, and1_out;
    two_in_and_gate and0(.a(a_old), .b(b_new), .z(and0_out));
    two_in_and_gate and1(.a(a_new), .b(b_old), .z(and1_out));
    tri_state tri_0(.a(and0_out), .enable(b_new), .out(FA_in0));
    tri_state tri_1(.a(and1_out), .enable(b_new), .out(FA_in1));
    tri_state tri_2(.a(carry_in), .enable(b_new), .out(FA_carry_in));
    Full_Adder FA(.a(FA_in0), .b(FA_in1), .carry_in(FA_carry_in), .sum(FA_out1), .carry_out(FA_out0));
    mux mux0(.a(last_carry), .b(FA_out0), .select(b_new), .out(carry_out));
    mux mux1(.a(and1_out), .b(FA_out1), .select(b_new), .out(sum));
endmodule

module Full_Adder_RB_three_In(input last_carry, a_old, b_new, in1, carry_in, output carry_out, sum);
    wire FA_in0, FA_in1, FA_carry_in, FA_out0, FA_out1, and0_out;
    two_in_and_gate and0(.a(a_old), .b(b_new), .z(and0_out));
    tri_state tri_0(.a(and0_out), .enable(b_new), .out(FA_in0));
    tri_state tri_1(.a(in1), .enable(b_new), .out(FA_in1));
    tri_state tri_2(.a(carry_in), .enable(b_new), .out(FA_carry_in));
    Full_Adder FA(.a(FA_in0), .b(FA_in1), .carry_in(FA_carry_in), .sum(FA_out1), .carry_out(FA_out0));
    mux mux0(.a(last_carry), .b(FA_out0), .select(b_new), .out(carry_out));
    mux mux1(.a(in1), .b(FA_out1), .select(b_new), .out(sum));
endmodule

module ArrayMultiplier_row_bypass_generic(product, a, b);
    parameter m = 64;
    parameter n = 64;
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
            //Full_Adder_RB_Four_In(input last_carry, a_old, a_new, b_old, b_new, carry_in, output carry_out, sum);
            Full_Adder_RB_Four_In FA_first(.last_carry(1'b0), .a_old(a[i]), .a_new(a[i+1]), .b_old(b[0]), .b_new(b[1]), .carry_in(1'b0),
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
                //Full_Adder_RB_three_In(input last_carry, a_old, b_new, in1, carry_in, output carry_out, sum);
                Full_Adder_RB_three_In FA_middle(.last_carry(c_partial[((m-1)*k)+j+1]), .a_old(a[j]), .b_new(b[k+2]), .in1(s_partial[((m-1)*k)+(j+1)]), .carry_in(c_partial[((m-1)*k)+j]),
                 .sum(s_partial[((m-1)*(k+1))+j]), .carry_out(c_partial[((m-1)*(k+1))+j]));
            end
        end
    endgenerate

    //generate middle lines of the multiplier - last column of each row
    genvar z;
    generate
        for(z=0; z<n-2; z=z+1)
        begin
            //Full_Adder_RB_Four_In(input last_carry, a_old, a_new, b_old, b_new, carry_in, output carry_out, sum);
            Full_Adder_RB_Four_In FA_middle_last_column(.last_carry(1'b0), .a_old(a[m-2]), .a_new(a[m-1]), .b_old(b[z+1]), .b_new(b[z+2]), .carry_in(c_partial[(m-1)*(z)+(m-2)]),
             .sum(s_partial[((m-1)*(z+1))+(m-2)]), .carry_out(c_partial[((m-1)*(z+1))+(m-2)]));
        end
    endgenerate
    
    wire extra_logic_wires[2+(3*(n-3)):0];
    genvar h;
    generate
        for(h = 0; h<n-1; h=h+1)
        begin
            if(h==0) //first column of additional logic
            begin
                two_in_and_gate and_first_column(.a(c_partial[((m-1)*h)]), .b(b[h+2]), .z(extra_logic_wires[0]));
            end else if (h==1)
            begin
                Full_Adder FA_second_row_extra(.a(s_partial[((m-1)*h)]), .b(1'b0), .carry_in(extra_logic_wires[0]), 
                 .sum(extra_logic_wires[1]), .carry_out(extra_logic_wires[2]));
                two_in_and_gate and_second_row_extra(.a(c_partial[((m-1)*h)]), .b(b[h+2]), .z(extra_logic_wires[3]));
            end else if (h < m-2) //all rows after first two except last row
            begin
                Full_Adder FA_second_row_extra(.a(s_partial[((m-1)*h)]), .b(extra_logic_wires[3*(h-1)]), .carry_in(extra_logic_wires[3*(h-1)-1]), 
                 .sum(extra_logic_wires[3*(h)-2]), .carry_out(extra_logic_wires[3*(h)-1]));
                two_in_and_gate and_second_row_extra(.a(c_partial[((m-1)*h)]), .b(b[h+2]), .z(extra_logic_wires[3*h]));
            end
        end
    endgenerate

    //generate last line of the multiplier
    genvar l;
    generate
        for(l=0; l<n; l=l+1)
        begin
            if(l==0) //first FA in last row
            begin
                Full_Adder FA_LR_first_column(.a(s_partial[((m-1)*(n-2))]), .b(extra_logic_wires[3*(m-3)]), .carry_in(extra_logic_wires[3*(m-3)-1]),
                 .sum(extra_logic_wires[1+(3*(n-3))]), .carry_out(extra_logic_wires[2+(3*(n-3))]));
            end else if (l < n-1)
            begin
                if(l==1)
                begin
                    Full_Adder FA_LR(.a(s_partial[((m-1)*(n-2))+l]), .b(c_partial[((m-1)*(n-2))+l-1]), .carry_in(extra_logic_wires[2+(3*(n-3))]),
                     .sum(s_partial[((m-1)*(n-1))+l-1]), .carry_out(c_partial[((m-1)*(n-1))+l-1]));
                end else begin
                    Full_Adder FA_LR(.a(s_partial[((m-1)*(n-2))+l]), .b(c_partial[((m-1)*(n-2))+l-1]), .carry_in(extra_logic_wires[2+(3*(n-3))]),
                     .sum(s_partial[((m-1)*(n-1))+l-1]), .carry_out(c_partial[((m-1)*(n-1))+l-1]));
                end
            end else begin
                Full_Adder_Last FA_Last(.a_row(a[m-1]), .b_row(b[n-1]), .b(c_partial[((m-1)*(n-2))+l-1]), .carry_in(c_partial[((m-1)*(n-1))+1]),
                 .sum(s_partial[((m-1)*(n-1))+l-1]), .carry_out(c_partial[((m-1)*(n-1))+l-1]));
            end
        end
    endgenerate

    //end product bits from first and middle cells
    generate
        for(i=0; i<n-1; i=i+1)
        begin
            if(i==0)
            begin
                buf (product[i+1], s_partial[(m-1)*i]);
            end else begin
                buf (product[i+1], extra_logic_wires[(i*3)-2]);
            end
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