`timescale 1ns / 1ps

module Full_Adder_four_in(input a_old, a_new, b_old, b_new, carry_in, output reg sum, carry_out);
    always @(a_old or a_new or b_old or b_new or carry_in)
    begin
        sum = (a_old & b_new) ^ (a_new & b_old) ^ carry_in;
        carry_out = (a_old & b_new) & (a_new & b_old) | ((a_old & b_new)^(a_new & b_old)) & carry_in;
    end
endmodule

module Full_Adder_three_in(input a_this, b_row, last_sum, carry_in, output reg sum, carry_out);
    always @(a_this or b_row or last_sum or carry_in)
    begin
        sum = (a_this & b_row) ^ last_sum ^ carry_in;
        carry_out = (a_this & b_row) & last_sum | ((a_this & b_row)^last_sum) & carry_in;
    end
endmodule

module Full_Adder_two_in(input a, b, carry_in, output reg sum, carry_out);
    always @(a or b or carry_in)
    begin
        sum = a ^ b ^ carry_in;
        carry_out = a & b | (a^b) & carry_in;
    end
endmodule

module two_in_and_gate(input a, b, output reg z);
    always @(a or b)
    begin
        z = a & b;
    end
endmodule

module ArrayMultiplier_generic(product, a, b);
    parameter m = 4;
    parameter n = 4;
    output wire [m+n-1:0] product;
    input[m-1:0] a;
    input[n-1:0] b;

    wire c_partial[(m-1)*(n):0];
    wire s_partial[(m-1)*(n):0];

    //generate first product bit
    two_in_and_gate two_and(.a(a[0]), .b(b[0]), .z(product[0]));

    //generate first row of the multiplier
    genvar i;
    generate
        for(i=0; i<m-1; i=i+1)
        begin
                Full_Adder_four_in FA_first(.a_old(a[i]), .a_new(a[i+1]), .b_old(b[0]), .b_new(b[1]), .carry_in(1'b0), .sum(s_partial[i]), .carry_out(c_partial[i]));
        end
    endgenerate

    //generate middle rows of the multiplier except last column
    genvar j, k;
    generate
        for(k=0; k<n-2; k=k+1) //row
        begin
            for(j=0; j<m-2; j=j+1)  //column
            begin
                Full_Adder_three_in FA_middle(.a_this(a[j]), .b_row(b[k+2]), .last_sum(s_partial[((m-1)*k)+(j+1)]), .carry_in(c_partial[((m-1)*k)+j]),
                 .sum(s_partial[((m-1)*(k+1))+j]), .carry_out(c_partial[((m-1)*(k+1))+j]));
            end
        end
    endgenerate

    //generate middle lines of the multiplier - last column of each row
    genvar z;
    generate
        for(z=0; z<n-2; z=z+1)
        begin
            Full_Adder_four_in FA_middle_last_column(.a_old(a[m-2]), .a_new(a[m-1]), .b_old(b[z+1]), .b_new(b[z+2]), .carry_in(c_partial[(m-1)*(z)+(m-2)]),
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
                Full_Adder_two_in Last_row_first_column(.a(s_partial[((m-1)*(n-2))+1]), .b(c_partial[((m-1)*(n-2))+l]), .carry_in(1'b0),
                 .sum(s_partial[((m-1)*(n-1))+l]), .carry_out(c_partial[((m-1)*(n-1))+l]));
            end else if (l < m-2)
            begin
                Full_Adder_two_in Last_row_middle_columns(.a(s_partial[((m-1)*(n-2))+l+1]), .b(c_partial[((m-1)*(n-2))+l]), .carry_in(c_partial[((m-1)*(n-1))+l-1]),
                 .sum(s_partial[((m-1)*(n-1))+l]), .carry_out(c_partial[((m-1)*(n-1))+l]));
            end else begin
                Full_Adder_three_in Last_row_last_column(.a_this(a[m-1]), .b_row(b[n-1]), .last_sum(c_partial[((m-1)*(n-2))+l]), .carry_in(c_partial[((m-1)*(n-1))+l-1]),
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