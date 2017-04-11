//Copyright 2017 Andrey S. Ionisyan (anserion@gmail.com)
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

unit Unit1;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, FileUtil, Forms, Controls, Graphics, Dialogs, Grids,
  StdCtrls;

type

  { TForm1 }

  TForm1 = class(TForm)
    BTN_exit: TButton;
    BTN_save: TButton;
    BTN_load: TButton;
    BTN_start: TButton;
    BTN_reset: TButton;
    CB_gates: TComboBox;
    EDIT_qgates_num: TEdit;
    EDIT_qbits_num: TEdit;
    EDIT_step: TEdit;
    EDIT_observe: TEdit;
    Label1: TLabel;
    Label2: TLabel;
    Label3: TLabel;
    Label4: TLabel;
    Label5: TLabel;
    OpenDialog1: TOpenDialog;
    SaveDialog1: TSaveDialog;
    SG_cur_state: TStringGrid;
    SG_all_states: TStringGrid;
    SG_alg: TStringGrid;
    procedure BTN_exitClick(Sender: TObject);
    procedure BTN_loadClick(Sender: TObject);
    procedure BTN_resetClick(Sender: TObject);
    procedure BTN_saveClick(Sender: TObject);
    procedure BTN_startClick(Sender: TObject);
    procedure CB_gatesChange(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure SG_algSelectCell(Sender: TObject; aCol, aRow: Integer;
      var CanSelect: Boolean);
    procedure SG_all_statesSelectCell(Sender: TObject; aCol, aRow: Integer;
      var CanSelect: Boolean);
  private
    { private declarations }
  public
    { public declarations }
    procedure Reset_emul;
    procedure Show_step(cur_step:integer);
    procedure Show_all_states;
  end;

var
  Form1: TForm1;

implementation

{$R *.lfm}

uses math;

const sqrt2=1.414213562;
      sqrt3=1.732050808;
      sqrt5=2.236067977;

type
TComplex=record
   re,im:real;
end;

TComplexVector=array of TComplex;
TComplexMatrix=array of TComplexVector;

TIntegerComplex=record
   re,im:integer;
end;
TIntegerComplexVector=array of TIntegerComplex;
TIntegerComplexMatrix=array of TIntegerComplexVector;

type
   tbit=(zero,one);
   tbit_vector=array of tbit;
   tbit_table=array of tbit_vector;
   tqbit=record p0,p1: tcomplex; end;
   tqregister=TComplexVector;

   //====================================================================
   //битовые операции
   function get_bit(n,k:integer):integer;
   begin get_bit:=(n shr k) and 1; end;

   function get_bits(n,k,nn:integer):integer;
   var mask:integer;
   begin mask:=(1 shl nn)-1; get_bits:=(n shr k) and mask; end;

   function set_bit(n,k,value:integer):integer;
   begin set_bit:=(n and (not (1 shl k))) or (value shl k); end;

   function insert_bit(n,k,value:integer):integer;
   var mask:integer;
   begin
      mask:=(1 shl k)-1;
      insert_bit:=((((n shr k)shl 1)or value) shl k) or (n and mask);
   end;

   function insert_bits(n,k,value,n_value:integer):integer;
   var mask:integer;
   begin
      mask:=(1 shl k)-1;
      insert_bits:=((((n shr k)shl n_value)or value) shl k) or (n and mask);
   end;

   function delete_bit(n,k:integer):integer;
   var mask:integer;
   begin
      mask:=(1 shl k)-1;
      delete_bit:=((n shr (k+1))shl k) or (n and mask);
   end;

   //функция вычисляет значение 2^n
   function pow2(n:integer):integer; begin pow2:=1 shl n; end;

   //функция вычисляет ближайшее меньшее целое значение log_2(n)
   function log2(n:integer):integer;
   var res:integer;
   begin
      res:=0;
      while n>1 do begin n:=n>>1; res:=res+1; end;
      log2:=res;
   end;

   //функция выполняет округление числа N вверх до ближайшей степени двойки
   function Power2RoundUp(N:integer):integer;
   var NN:integer;
   begin
     NN:=(1<<log2(N)); if N>NN then NN:=NN<<1;
     Power2RoundUp:=NN;
   end;

   //====================================================================
function c_complex(re,im:real):TComplex;
begin c_complex.re:=re; c_complex.im:=im; end;

function c_vector(n:integer):TComplexVector;
var tmp:TComplexVector; begin setlength(tmp,n); c_vector:=tmp; end;

function c_matrix(n,m:integer):TComplexMatrix;
var i:integer; tmp:TComplexMatrix;
begin
   setlength(tmp,n);
   for i:=0 to n-1 do setlength(tmp[i],m);
   c_matrix:=tmp;
end;

procedure c_matrix_destroy(var M:TComplexMatrix);
var i,n:integer;
begin
   n:=length(M); for i:=0 to n-1 do setlength(M[i],0);
   setlength(M,0);
end;

function c_zero:TComplex; begin c_zero.re:=0; c_zero.im:=0; end;
function c_one:TComplex; begin c_one.re:=1; c_one.im:=0; end;
function c_minus_one:TComplex; begin c_minus_one.re:=-1; c_minus_one.im:=0; end;
function c_i:TComplex; begin c_i.re:=0; c_i.im:=1; end;
function c_minus_i:TComplex; begin c_minus_i.re:=0; c_minus_i.im:=-1; end;
function c_pi:TComplex; begin c_pi.re:=Pi; c_pi.im:=0; end;
function c_pi2:TComplex; begin c_pi2.re:=0.5*Pi; c_pi2.im:=0; end;
function c_pi4:TComplex; begin c_pi4.re:=0.25*Pi; c_pi4.im:=0; end;
function c_pi8:TComplex; begin c_pi8.re:=0.125*Pi; c_pi8.im:=0; end;
function c_sqrt2:TComplex; begin c_sqrt2.re:=sqrt2; c_sqrt2.im:=0; end;
function c_sqrt3:TComplex; begin c_sqrt3.re:=sqrt3; c_sqrt3.im:=0; end;
function c_sqrt5:TComplex; begin c_sqrt5.re:=sqrt5; c_sqrt5.im:=0; end;

function c_root_of_one_CCW(k,n:integer):TComplex;
var phi:real;
begin
     phi:=2*PI*k/n;
     c_root_of_one_CCW.re:=cos(phi);
     c_root_of_one_CCW.im:=sin(phi);
end;

function c_root_of_one_CW(k,n:integer):TComplex;
var phi:real;
begin
     phi:=-2*PI*k/n;
     c_root_of_one_CW.re:=cos(phi);
     c_root_of_one_CW.im:=sin(phi);
end;
//====================================================================

function c_dup(value:TComplex):TComplex;
begin c_dup.re:=value.re; c_dup.im:=value.im; end;

function c_conj(value:TComplex):TComplex;
begin c_conj.re:=value.re; c_conj.im:=-value.im; end;

function c_amp2(a:TComplex):real;
begin c_amp2:=sqr(a.re)+sqr(a.im); end;

function c_amp(a:TComplex):real;
begin c_amp:=sqrt(sqr(a.re)+sqr(a.im)); end;

function c_phi(a:TComplex):real;
var res,amp:real;
begin
     amp:=c_amp(a);
     if amp=0 then res:=0 else res:=arccos(a.re/amp);
     c_phi:=res;
end;

function c_TrigToAlg(amp,phi:real):TComplex;
begin c_TrigToAlg.re:=amp*cos(phi); c_TrigToAlg.im:=amp*sin(phi); end;

function c_amp_cmp(a,b:TComplex):integer;
var amp2_a,amp2_b:real; res:integer;
begin
     res:=0;
     amp2_a:=sqr(a.re)+sqr(a.im);
     amp2_b:=sqr(b.re)+sqr(b.im);
     if amp2_a>amp2_b then res:=1;
     if amp2_a=amp2_b then res:=0;
     if amp2_a<amp2_b then res:=-1;
     c_amp_cmp:=res;
end;

function c_neg(a:TComplex):TComplex;
begin c_neg.re:=-a.re; c_neg.im:=-a.im; end;

function c_inv(a:TComplex):TComplex;
var amp,inv_phi:real; res:TComplex;
begin
   amp:=c_amp(a);
   if amp=0 then res:=c_zero
      else
      begin
         inv_phi:=-arccos(a.re/amp);
         res:=c_TrigToAlg(1.0/amp,inv_phi);
      end;
   c_inv:=res;
end;

function c_add(a,b:TComplex):TComplex;
begin c_add.re:=a.re+b.re; c_add.im:=a.im+b.im; end;

function c_sub(a,b:TComplex):TComplex;
begin c_sub.re:=a.re-b.re; c_sub.im:=a.im-b.im; end;

function c_mul(a,b:TComplex):TComplex;
begin c_mul.re:=a.re*b.re-a.im*b.im; c_mul.im:=a.re*b.im+a.im*b.re; end;

function c_div(a,b:TComplex):TComplex;
begin
     c_div.re:=(a.re*b.re+a.im*b.im)/(b.re*b.re+b.im*b.im);
     c_div.im:=(a.im*b.re-a.re*b.im)/(b.re*b.re+b.im*b.im);
end;

procedure c_AlgToTrig(alg:TComplex; var amp,phi:real);
begin
     amp:=c_amp(alg);
     phi:=c_phi(alg);
end;

function c_sqr(arg:TComplex):TComplex;
begin c_sqr.re:=arg.re*arg.re-arg.im*arg.im; c_sqr.im:=2*arg.re*arg.im; end;

function c_exp_ix(x:real):TComplex;
begin c_exp_ix.re:=cos(x); c_exp_ix.im:=sin(x); end;

procedure c_sqrt(arg:TComplex; var res1,res2:TComplex);
var amp,phi:real;
begin
     amp:=sqrt(c_amp(arg));
     phi:=c_phi(arg);
     res1.re:=amp*cos(phi/2); res1.im:=amp*sin(phi/2);
     res2.re:=amp*cos(phi/2+PI); res2.im:=amp*sin(phi/2+PI);
end;

function c_exp(arg:TComplex):TComplex;
var exp_x:real;
begin
     exp_x:=exp(arg.re);
     c_exp.re:=exp_x*cos(arg.im);
     c_exp.im:=exp_x*sin(arg.im);
end;

function c_ln(arg:TComplex; k:integer):TComplex;
var amp,phi:real;
begin
     amp:=c_amp(arg);
     phi:=c_phi(arg);
     if amp>0 then amp:=ln(amp);
     c_ln.re:=amp;
     c_ln.im:=phi+2*PI*k;
end;

function c_power(arg,pow:TComplex; k:integer):TComplex;
begin
     c_power:=c_exp(c_mul(pow,c_ln(arg,k)));
end;

//-------------------------------------------------------------
procedure c_vector_fill(value:TComplex; var V:TComplexVector);
var i,n:integer;
begin n:=length(V); for i:=0 to n-1 do V[i]:=value; end;

procedure c_vector_copy(var src,dst:TComplexVector);
var i,n:integer;
begin n:=length(src); for i:=0 to n-1 do dst[i]:=src[i]; end;

procedure c_vectors_swap(var V1,V2:TComplexVector);
var i,n:integer; tmp:TComplex;
begin
   n:=length(V1);
   for i:=0 to n-1 do begin tmp:=V1[i]; V1[i]:=V2[i]; V2[i]:=tmp; end;
end;

procedure c_matrix_fill(value:TComplex; var A:TComplexMatrix);
var i,n:integer;
begin
   n:=length(A);
   for i:=0 to n-1 do c_vector_fill(value,A[i]);
end;

procedure c_matrix_copy(var src,dst:TComplexMatrix);
var i,n:integer;
begin n:=length(src); for i:=0 to n-1 do c_vector_copy(src[i],dst[i]); end;

procedure c_subvector_to_vector_put(k:integer; var subvector,V:TComplexVector);
var i,sn,n:integer;
begin
   sn:=length(subvector); n:=length(V);
   if k+sn>n then sn:=n-k;
   for i:=0 to sn-1 do V[i+k]:=subvector[i];
end;

procedure c_subvector_from_vector_get(k,n:integer; var subvector,V:TComplexVector);
var i,nn:integer;
begin
   nn:=length(V);
   if k+n>nn then n:=nn-n;
   for i:=0 to n-1 do subvector[i]:=V[k+i];
end;

procedure c_matrix_raw_put(k:integer; var A:TComplexMatrix; var V:TComplexVector);
var i,n:integer;
begin n:=length(V); for i:=0 to n-1 do A[k,i]:=V[i]; end;

procedure c_matrix_col_put(k:integer; var A:TComplexMatrix; var V:TComplexVector);
var i,m:integer;
begin m:=length(V); for i:=0 to m-1 do A[i,k]:=V[i]; end;

procedure c_matrix_raw_get(k:integer; var A:TComplexMatrix; var V:TComplexVector);
var i,n:integer;
begin n:=length(A[k]); for i:=0 to n-1 do V[i]:=A[k,i]; end;

procedure c_matrix_col_get(k:integer; var A:TComplexMatrix; var V:TComplexVector);
var i,m:integer;
begin m:=length(A); for i:=0 to m-1 do V[i]:=A[i,k]; end;

procedure c_submatrix_to_matrix_put(raw,col:integer; var submatrix,A:TComplexMatrix);
var i,j,n,m:integer;
begin
   n:=length(submatrix);  m:=length(submatrix[0]);
   for i:=0 to n-1 do
      for j:=0 to m-1 do
         A[raw+i,col+j]:=submatrix[i,j];
end;

procedure c_submatrix_from_matrix_get(raw,col,n,m:integer; var submatrix,A:TComplexMatrix);
var i,j:integer;
begin
   for i:=0 to n-1 do
      for j:=0 to m-1 do
         submatrix[i,j]:=A[raw+i,col+j];
end;

procedure c_matrix_transp(var A,res:TComplexMatrix);
var i,n,m:integer; tmp:TComplexVector;
begin
   n:=length(A); m:=length(A[0]); setlength(tmp,m);
   for i:=0 to n-1 do
   begin
      c_matrix_raw_get(i,A,tmp);
      c_matrix_col_put(i,res,tmp);
   end;
end;

procedure c_marix_hermitian(var A,res:TComplexMatrix);
var i,j,n,m:integer;
begin
   c_matrix_transp(A,res);
   n:=length(res); m:=length(res[0]);
   for i:=0 to n-1 do
   for j:=0 to m-1 do
      res[i,j].im:=-res[i,j].im;
end;
//-------------------------------------------------------------

function c_vector_summ(var a:TComplexVector):TComplex;
var tmp:TComplex; i,n:integer;
begin
   n:=length(a); tmp:=c_zero;
   for i:=0 to n-1 do tmp:=c_add(a[i],tmp);
   c_vector_summ:=tmp;
end;

function c_vector_prod(var a:TComplexVector):TComplex;
var tmp:TComplex; i,n:integer;
begin
   n:=length(a); tmp:=c_one;
   for i:=0 to n-1 do tmp:=c_mul(a[i],tmp);
   c_vector_prod:=tmp;
end;

function c_vector_dist2(var a:TComplexVector):TComplex;
var tmp:TComplex; i,n:integer;
begin
   n:=length(a); tmp:=c_zero;
   for i:=0 to n-1 do tmp:=c_add(c_mul(a[i],a[i]),tmp);
   c_vector_dist2:=tmp;
end;

function c_vector_mean(var a:TComplexVector):TComplex;
begin c_vector_mean:=c_div(c_vector_summ(a),c_complex(length(a),0)); end;

function c_vector_dispersion(var a:TComplexVector):TComplex;
var i,n:integer; average,tmp:TComplex;
begin
   n:=length(a);
   average:=c_vector_mean(a);
   tmp:=c_zero;
   for i:=0 to n-1 do tmp:=c_add(tmp,c_sqr(c_sub(a[i],average)));
   c_vector_dispersion:=c_div(tmp,c_complex(n,0));
end;

procedure c_vector_diff(var a,res:TComplexVector);
var i,n:integer; c_two:TComplex;
begin
   n:=length(a); c_two:=c_complex(2,0);
   res[0]:=c_sub(a[1],a[0]); res[n-1]:=c_sub(a[n-1],a[n-2]);
   for i:=1 to n-2 do
      res[i]:=c_div(c_sub(a[i+1],a[i-1]),c_two);
end;

procedure c_func_diff(var x,f,res:TComplexVector);
var i,n:integer;
begin
   n:=length(x);
   res[0]:=c_div(c_sub(f[1],f[0]),c_sub(x[1],x[0]));
   res[n-1]:=c_div(c_sub(f[n-1],f[n-2]),c_sub(x[n-1],x[n-2]));
   for i:=1 to n-2 do
      res[i]:=c_div(c_sub(f[i+1],f[i-1]),c_sub(x[i+1],x[i-1]));
end;

procedure c_vector_neg(var a,neg_a:TComplexVector);
var i,n:integer;
begin n:=length(a); for i:=0 to n-1 do neg_a[i]:=c_neg(a[i]); end;

procedure c_vector_add_vector(var a,b,c:TComplexVector);
var i,n:integer;
begin n:=length(a); for i:=0 to n-1 do c[i]:=c_add(a[i],b[i]); end;

procedure c_vector_sub_vector(var a,b,c:TComplexVector);
var i,n:integer;
begin n:=length(a); for i:=0 to n-1 do c[i]:=c_sub(a[i],b[i]); end;

function c_vectors_scalar_mul(var a,b:TComplexVector):TComplex;
var tmp:TComplex; i,n:integer;
begin
   n:=length(a); tmp:=c_zero;
   for i:=0 to n-1 do tmp:=c_add(c_mul(a[i],b[i]),tmp);
   c_vectors_scalar_mul:=tmp;
end;

function c_vectors_convolution(var a,b:TComplexVector):TComplex;
var tmp:TComplex; i,n:integer;
begin
   n:=length(a); tmp:=c_zero;
   for i:=0 to n-1 do tmp:=c_add(c_mul(a[i],b[n-i-1]),tmp);
   c_vectors_convolution:=tmp;
end;

procedure c_vector_add_scalar(lambda:TComplex; var a,res:TComplexVector);
var i,n:integer;
begin n:=length(a); for i:=0 to n-1 do res[i]:=c_add(lambda,a[i]); end;

procedure c_vector_mul_scalar(lambda:TComplex; var a,res:TComplexVector);
var i,n:integer;
begin n:=length(a); for i:=0 to n-1 do res[i]:=c_mul(lambda,a[i]); end;

procedure c_matrix_mul_scalar(lambda:TComplex; var a,res:TComplexMatrix);
var i,j,n,m:integer;
begin
   n:=length(a); m:=length(a[0]);
   for i:=0 to n-1 do
      for j:=0 to m-1 do
         res[i,j]:=c_mul(lambda,a[i,j]);
end;

procedure c_matrix_add_matrix(var a,b,c:TComplexMatrix);
var i,j,n,m:integer;
begin
   n:=length(a); m:=length(a[0]);
   for i:=0 to n-1 do
      for j:=0 to m-1 do
         c[i,j]:=c_add(a[i,j],b[i,j]);
end;

procedure c_matrix_sub_matrix(var a,b,c:TComplexMatrix);
var i,j,n,m:integer;
begin
   n:=length(a); m:=length(a[0]);
   for i:=0 to n-1 do
      for j:=0 to m-1 do
         c[i,j]:=c_sub(a[i,j],b[i,j]);
end;

procedure c_matrix_mul_vector(var A:TComplexMatrix; var V,res:TComplexVector);
var i,n:integer;
begin
   n:=length(A);
   for i:=0 to n-1 do res[i]:=c_vectors_scalar_mul(A[i],V);
end;

procedure c_vector_mul_matrix(var V:TComplexVector; var A:TComplexMatrix; var res:TComplexVector);
var i,n,m:integer; tmp:TComplexVector;
begin
   m:=length(V); n:=length(A); setlength(tmp,n);
   for i:=0 to m-1 do
   begin
      c_matrix_col_get(i,A,tmp);
      res[i]:=c_vectors_scalar_mul(V,tmp);
   end;
end;

procedure c_matrix_mul_matrix(var A,B,C:TComplexMatrix);
var i,n:integer;
begin
   n:=length(A);
   for i:=0 to n-1 do c_vector_mul_matrix(A[i],B,C[i]);
end;

procedure c_vectorv_kronmul_vectorh(var a,b:TComplexVector; var c:TComplexMatrix);
var i,j,n,m:integer;
begin
   n:=length(a); m:=length(b);
   for i:=0 to n-1 do
      for j:=0 to m-1 do
         c[i,j]:=c_mul(a[i],b[j]);
end;

procedure c_vectorh_kronmul_vectorv(var a,b:TComplexVector; var c:TComplexMatrix);
var i,j,n,m:integer;
begin
   n:=length(b); m:=length(a);
   for i:=0 to n-1 do
      for j:=0 to m-1 do
         c[i,j]:=c_mul(a[j],b[i]);
end;

procedure c_matrix_kronmul_matrix(var a,b,c:TComplexMatrix);
var ia,ja,ib,jb,na,ma,nb,mb:integer;
begin
   na:=length(a); ma:=length(a[0]);
   nb:=length(b); mb:=length(b[0]);
   for ia:=0 to na-1 do
   for ja:=0 to ma-1 do
      for ib:=0 to nb-1 do
      for jb:=0 to mb-1 do
         c[ia*nb+ib,ja*mb+jb]:=c_mul(a[ia,ja],b[ib,jb]);
end;

//====================================================================

function qregister_create(n:integer):tqregister;
var tmp:tqregister;
begin
   setlength(tmp,pow2(n));
   c_vector_fill(c_zero,tmp);
   qregister_create:=tmp;
end;

procedure qregister_set_tbits(var qreg:tqregister; value:tbit_vector);
var i,n,idx:integer;
begin
   c_vector_fill(c_zero,qreg);
   n:=length(value);
   if value[n-1]=one then idx:=1 else idx:=0;
   for i:=n-2 downto 0 do
   begin
      idx:=idx*2;
      if value[i]=one then idx:=idx+1;
   end;
   qreg[idx]:=c_one;
end;

procedure qregister_get_tbits(var qreg:tqregister; res:tbit_vector);
var i,n,nn,max_amp2_idx:integer; amp2,max_amp2:real;
begin
   nn:=length(qreg);
   max_amp2_idx:=0; max_amp2:=c_amp2(qreg[0]);
   for i:=1 to nn-1 do
   begin
      amp2:=c_amp2(qreg[i]);
      if amp2>max_amp2 then begin max_amp2:=amp2; max_amp2_idx:=i; end;
   end;
   n:=log2(nn);
   for i:=0 to n-1 do
   begin
      if (max_amp2_idx mod 2)=0 then res[i]:=zero else res[i]:=one;
      max_amp2_idx:=max_amp2_idx div 2;
   end;
end;

procedure qregister_apply_gate(var qreg:tqregister; gate:TComplexMatrix);
var tmp_reg:tqregister;
begin
   tmp_reg:=c_vector(length(qreg));
   c_vector_copy(qreg,tmp_reg);
   c_matrix_mul_vector(gate,qreg,tmp_reg);
   c_vector_copy(tmp_reg,qreg);
   setlength(tmp_reg,0);
end;

function qregister_calc_gate1(k,n:integer; p00,p01,p10,p11:tcomplex):TComplexMatrix;
var i:integer;
   tmp,tmp1,gateP,gateI:TComplexMatrix;
begin
   gateI:=c_matrix(2,2);
   gateI[0,0]:=c_one; gateI[0,1]:=c_zero; gateI[1,0]:=c_zero; gateI[1,1]:=c_one;
   gateP:=c_matrix(2,2);
   gateP[0,0]:=p00; gateP[0,1]:=p01; gateP[1,0]:=p10; gateP[1,1]:=p11;
   tmp:=c_matrix(2,2);
   if k=0 then c_matrix_copy(gateP,tmp) else c_matrix_copy(gateI,tmp);
   for i:=1 to n-1 do
   begin
      c_matrix_destroy(tmp1);
      tmp1:=c_matrix(pow2(i+1),pow2(i+1));
      if i=k then c_matrix_kronmul_matrix(gateP,tmp,tmp1)
             else c_matrix_kronmul_matrix(gateI,tmp,tmp1);
      c_matrix_destroy(tmp);
      tmp:=c_matrix(pow2(i+1),pow2(i+1));
      c_matrix_copy(tmp1,tmp);
   end;
   c_matrix_destroy(tmp1);
   c_matrix_destroy(gateI);
   c_matrix_destroy(gateP);
   qregister_calc_gate1:=tmp;
end;

function qregister_calc_I_gate(n:integer):TComplexMatrix;
var i,nn:integer; tmp:TComplexMatrix;
begin
   nn:=pow2(n);
   tmp:=c_matrix(nn,nn);
   for i:=0 to nn-1 do tmp[i,i]:=c_one;
   qregister_calc_I_gate:=tmp;
end;

function qregister_calc_H_gate(n:integer):TComplexMatrix;
var i:integer; tmp,tmp1,gateH:TComplexMatrix;
begin
   gateH:=c_matrix(2,2);
   gateH[0,0]:=c_complex(1.0/sqrt2,0);
   gateH[0,1]:=c_complex(1.0/sqrt2,0);
   gateH[1,0]:=c_complex(1.0/sqrt2,0);
   gateH[1,1]:=c_complex(-1.0/sqrt2,0);
   tmp:=c_matrix(2,2); c_matrix_copy(gateH,tmp);
   for i:=1 to n-1 do
   begin
      c_matrix_destroy(tmp1);
      tmp1:=c_matrix(pow2(i+1),pow2(i+1));
      c_matrix_kronmul_matrix(gateH,tmp,tmp1);
      c_matrix_destroy(tmp);
      tmp:=c_matrix(pow2(i+1),pow2(i+1));
      c_matrix_copy(tmp1,tmp);
   end;
   c_matrix_destroy(tmp1);
   c_matrix_destroy(gateH);
   qregister_calc_H_gate:=tmp;
end;

function qregister_calc_X_gate(k,n:integer):TComplexMatrix;
begin
qregister_calc_X_gate:=qregister_calc_gate1(k,n,c_zero,c_one,c_one,c_zero);
end;

function qregister_calc_Y_gate(k,n:integer):TComplexMatrix;
begin
qregister_calc_Y_gate:=qregister_calc_gate1(k,n,c_zero,c_minus_i,c_i,c_zero);
end;

function qregister_calc_Z_gate(k,n:integer):TComplexMatrix;
begin
qregister_calc_Z_gate:=qregister_calc_gate1(k,n,c_one,c_zero,c_zero,c_minus_one);
end;

function qregister_calc_PHI_gate(k,n:integer; phi:real):TComplexMatrix;
begin
qregister_calc_PHI_gate:=qregister_calc_gate1(k,n,c_one,c_zero,c_zero,c_exp_ix(phi));
end;

function qregister_calc_NOT_gate(k,n:integer):TComplexMatrix;
var i,nn,nn2,idx1,idx2:integer; tmp:TComplexMatrix;
begin
   nn:=pow2(n); nn2:=nn shr 1;
   tmp:=qregister_calc_I_gate(n);
   for i:=0 to nn2-1 do
   begin
      idx1:=insert_bit(i,k,0);
      idx2:=insert_bit(i,k,1);
      c_vectors_swap(tmp[idx1],tmp[idx2]);
   end;
   qregister_calc_NOT_gate:=tmp;
end;

function qregister_calc_SWAP_gate(c1_qbit,c2_qbit,n:integer):TComplexMatrix;
var i,nn,idx,b1,b2:integer; tmp:TComplexMatrix;
begin
   nn:=pow2(n);
   tmp:=c_matrix(nn,nn); c_matrix_fill(c_zero,tmp);
   for i:=0 to nn-1 do
   begin
      b1:=get_bit(i,c1_qbit); b2:=get_bit(i,c2_qbit);
      idx:=set_bit(i,c1_qbit,b2); idx:=set_bit(idx,c2_qbit,b1);
      tmp[i,idx]:=c_one;
   end;
   qregister_calc_SWAP_gate:=tmp;
end;

function qregister_calc_CNOT_gate(c_qbit,u_qbit,n:integer):TComplexMatrix;
var i,nn,nn2,idx1,idx2:integer; tmp:TComplexMatrix;
begin
   nn:=pow2(n); nn2:=nn shr 1;
   tmp:=qregister_calc_I_gate(n);
   for i:=0 to nn2-1 do
   begin
      idx1:=insert_bit(i,u_qbit,0);
      idx2:=insert_bit(i,u_qbit,1);
      if get_bit(idx1,c_qbit)=1 then c_vectors_swap(tmp[idx1],tmp[idx2]);
   end;
   qregister_calc_CNOT_gate:=tmp;
end;

procedure qregister_get_qbits(start_pos,n_res:integer; var qreg,res:tqregister);
var n,i,idx:integer;
begin
   n:=length(qreg);
   c_vector_fill(c_zero,res);
   for i:=0 to n-1 do
   begin
      idx:=get_bits(i,start_pos,n_res);
      res[idx]:=c_add(res[idx],qreg[i]);
   end;
end;

procedure qregister_erase_qbits(start_pos,n_erase:integer; var qreg,res:tqregister);
var n,nq,i,idx,end_pos,left_part,right_part:integer;
begin
   n:=length(qreg); nq:=log2(n);
   end_pos:=start_pos+n_erase;
   c_vector_fill(c_zero,res);
   for i:=0 to n-1 do
   begin
      left_part:=get_bits(i,end_pos,nq-end_pos);
      right_part:=get_bits(i,0,start_pos);
      idx:=insert_bits(left_part,0,right_part,start_pos);
      res[idx]:=c_add(res[idx],qreg[i]);
   end;
end;

procedure qregister_increase_qbits_num(n:integer; var qreg,res:tqregister);
var i,nn_qreg,nn_res:integer;
begin
   nn_qreg:=length(qreg);
   nn_res:=pow2(log2(nn_qreg)+n);
   for i:=0 to nn_qreg-1 do res[i]:=qreg[i];
   for i:=nn_qreg to nn_res-1 do res[i]:=c_zero;
end;

//====================================================================

function dec_to_bin(value:integer):string;
var s:string;
begin
     if value=0 then s:='0' else s:='';
     while value>0 do
     begin
        if (value and 1)=1 then s:='1'+s else s:='0'+s;
        value:=value shr 1;
     end;
     dec_to_bin:=s;
end;

function dec_to_nbin(value,n:integer):string;
var i:integer; s:string;
begin
     s:='';
     for i:=1 to n do
     begin
        if (value and 1)=1 then s:='1'+s else s:='0'+s;
        value:=value shr 1;
     end;
     dec_to_nbin:=s;
end;

//======================================================================
var
   q_gates:array of TComplexMatrix;
   q_regs:array of tqregister;
   qbits_num,qgates_num,states_num:integer;
//======================================================================

procedure TForm1.Show_step(cur_step:integer);
var
   state_idx,qbit_idx:integer;
   observe_bits:tbit_vector;
   s:string;
begin
  EDIT_step.text:=IntToStr(cur_step);

  SG_cur_state.Cells[0,0]:='N';
  SG_cur_state.Cells[1,0]:='bin';
  SG_cur_state.Cells[2,0]:='RE';
  SG_cur_state.Cells[3,0]:='IM';
  SG_cur_state.Cells[4,0]:='AMP';
  SG_cur_state.Cells[5,0]:='PHI';
  for state_idx:=1 to states_num do
  begin
    SG_cur_state.Cells[0,state_idx]:=IntToStr(state_idx-1);
    SG_cur_state.Cells[1,state_idx]:=Dec_to_nbin(state_idx-1,qbits_num);
    SG_cur_state.Cells[2,state_idx]:=FloatToStr(q_regs[cur_step,state_idx-1].re);
    SG_cur_state.Cells[3,state_idx]:=FloatToStr(q_regs[cur_step,state_idx-1].im);
    SG_cur_state.Cells[4,state_idx]:=FloatToStr(c_amp(q_regs[cur_step,state_idx-1]));
    SG_cur_state.Cells[5,state_idx]:=FloatToStr(c_phi(q_regs[cur_step,state_idx-1])*180/Pi);
  end;

  setlength(observe_bits,qbits_num);
  qregister_get_tbits(q_regs[cur_step],observe_bits);
  s:='';
  for qbit_idx:=0 to qbits_num-1 do
      if observe_bits[qbit_idx]=zero then s:='0'+s else s:='1'+s;
  EDIT_observe.text:=s;
  setlength(observe_bits,0);
end;

procedure TForm1.Show_all_states;
var gate_idx,state_idx:integer;
begin
  SG_all_states.Cells[0,0]:='N';
  for state_idx:=1 to states_num do SG_all_states.Cells[0,state_idx]:=IntToStr(state_idx-1);
  for gate_idx:=1 to qgates_num do SG_all_states.Cells[gate_idx,0]:=IntToStr(gate_idx);
  for state_idx:=1 to states_num do
      for gate_idx:=1 to qgates_num do
          SG_all_states.Cells[gate_idx,state_idx]:=FloatToStr(c_amp(q_regs[gate_idx,state_idx-1]));
end;

procedure TForm1.reset_emul;
var gate_idx,qbit_idx:integer;
begin
  for gate_idx:=1 to length(q_gates)-1 do
  begin
    if length(q_gates[gate_idx])<>0 then c_matrix_destroy(q_gates[gate_idx]);
    if length(q_regs[gate_idx])<>0 then setlength(q_regs[gate_idx],0);
  end;

  qgates_num:=StrToInt(EDIT_qgates_num.text);
  qbits_num:=StrToInt(EDIT_qbits_num.text);
  states_num:=pow2(qbits_num);

  setlength(q_gates,qgates_num+1);
  setlength(q_regs,qgates_num+1);

  for gate_idx:=0 to qgates_num do
  begin
    q_gates[gate_idx]:=qregister_calc_I_gate(qbits_num);
    q_regs[gate_idx]:=qregister_create(qbits_num);
  end;

  SG_cur_state.RowCount:=states_num+1;

  SG_all_states.RowCount:=states_num+1;
  SG_all_states.ColCount:=qgates_num+1;

  SG_alg.RowCount:=qbits_num+1;
  SG_alg.ColCount:=qgates_num+1;

  for qbit_idx:=1 to qbits_num do SG_alg.Cells[0,qbit_idx]:='qbit_'+IntToStr(qbit_idx-1);
  for gate_idx:=1 to qgates_num do SG_alg.Cells[gate_idx,0]:=IntToStr(gate_idx);
  for qbit_idx:=1 to qbits_num do SG_alg.Cells[1,qbit_idx]:='RESET';
  for gate_idx:=2 to qgates_num do
      for qbit_idx:=1 to qbits_num do SG_alg.Cells[gate_idx,qbit_idx]:='';

  Show_all_states;
  Show_step(1);

  CB_gates.Items.Clear;
  CB_gates.Height:=SG_alg.DefaultRowHeight;
  CB_gates.Items.Add('');
  CB_gates.Items.Add('SET');
  CB_gates.Items.Add('RESET');
  CB_gates.Items.Add('I');
  CB_gates.Items.Add('H');
  CB_gates.Items.Add('X');
  CB_gates.Items.Add('Y');
  CB_gates.Items.Add('Z');
  CB_gates.Items.Add('NOT');
  CB_gates.Items.Add('PHI_p30');
  CB_gates.Items.Add('PHI_p45');
  CB_gates.Items.Add('PHI_p60');
  CB_gates.Items.Add('PHI_p90');
  CB_gates.Items.Add('PHI_180');
  CB_gates.Items.Add('PHI_m30');
  CB_gates.Items.Add('PHI_m45');
  CB_gates.Items.Add('PHI_m60');
  CB_gates.Items.Add('PHI_m90');
  CB_gates.Items.Add('CNOT_C');
  CB_gates.Items.Add('CNOT_U');
  CB_gates.Items.Add('SWAP1');
  CB_gates.Items.Add('SWAP2');
  CB_gates.ItemIndex:=0;
  CB_gates.visible:=false;
end;

{ TForm1 }

procedure TForm1.BTN_exitClick(Sender: TObject);
begin
  close;
end;

procedure TForm1.BTN_loadClick(Sender: TObject);
var f:TextFile; gate_idx,qbit_idx,char_idx:integer; gate_name,s:string;
begin
  if OpenDialog1.Execute then
  begin
     AssignFile(f,OpenDialog1.filename);
     Reset(f);
     readln(f,qbits_num); EDIT_qbits_num.text:=IntToStr(qbits_num);
     readln(f,qgates_num); EDIT_qgates_num.text:=IntToStr(qgates_num);
     RESET_emul;
     while not(eof(f)) do
     begin
        read(f,gate_idx);
        if gate_idx<1 then gate_idx:=1;
        if gate_idx>qgates_num then gate_idx:=qgates_num;
        read(f,qbit_idx);
        if qbit_idx<0 then qbit_idx:=0;
        if qbit_idx>qbits_num-1 then qbit_idx:=qbits_num-1;
        readln(f,s);
        gate_name:='';
        for char_idx:=1 to length(s) do
            if s[char_idx]<>' ' then gate_name:=gate_name+s[char_idx];
        SG_alg.Cells[gate_idx,qbit_idx+1]:=gate_name;
     end;
     CloseFile(f);
  end;
end;

procedure TForm1.BTN_resetClick(Sender: TObject);
begin
  RESET_emul;
end;

procedure TForm1.BTN_saveClick(Sender: TObject);
var f:TextFile; gate_idx,qbit_idx:integer;
begin
  if SaveDialog1.Execute then
  begin
    AssignFile(f,SaveDialog1.filename);
    Rewrite(f);
    writeln(f,qbits_num);
    writeln(f,qgates_num);
    for gate_idx:=1 to qgates_num do
        for qbit_idx:=1 to qbits_num do
            if SG_alg.Cells[gate_idx,qbit_idx]<>'' then
               writeln(f,gate_idx,' ',qbit_idx-1,' ',SG_alg.Cells[gate_idx,qbit_idx]);
    CloseFile(f);
  end;
end;

procedure TForm1.BTN_startClick(Sender: TObject);
var
   bin_reg:tbit_vector;
   gate_idx,qbit_idx,state_idx:integer;
   swap1_idx,swap2_idx,cnot_c_idx,cnot_u_idx:integer;
begin
  setlength(bin_reg,qbits_num);
  for qbit_idx:=1 to qbits_num do
      if SG_alg.Cells[1,qbit_idx]='SET' then bin_reg[qbit_idx-1]:=one
                                        else bin_reg[qbit_idx-1]:=zero;
  qregister_set_tbits(q_regs[1],bin_reg);

  for gate_idx:=1 to qgates_num do c_matrix_destroy(q_gates[gate_idx]);

  for gate_idx:=2 to qgates_num do
  begin
      swap1_idx:=-1; swap2_idx:=-1;
      cnot_c_idx:=-1; cnot_u_idx:=-1;
      for qbit_idx:=1 to qbits_num do
      begin
         if SG_alg.Cells[gate_idx,qbit_idx]='I' then q_gates[gate_idx]:=qregister_calc_I_gate(qbits_num);
         if SG_alg.Cells[gate_idx,qbit_idx]='H' then q_gates[gate_idx]:=qregister_calc_H_gate(qbits_num);
         if SG_alg.Cells[gate_idx,qbit_idx]='X' then q_gates[gate_idx]:=qregister_calc_X_gate(qbit_idx-1,qbits_num);
         if SG_alg.Cells[gate_idx,qbit_idx]='Y' then q_gates[gate_idx]:=qregister_calc_Y_gate(qbit_idx-1,qbits_num);
         if SG_alg.Cells[gate_idx,qbit_idx]='Z' then q_gates[gate_idx]:=qregister_calc_Z_gate(qbit_idx-1,qbits_num);
         if SG_alg.Cells[gate_idx,qbit_idx]='NOT' then q_gates[gate_idx]:=qregister_calc_NOT_gate(qbit_idx-1,qbits_num);
         if SG_alg.Cells[gate_idx,qbit_idx]='PHI_p30' then q_gates[gate_idx]:=qregister_calc_PHI_gate(qbit_idx-1,qbits_num,PI/6);
         if SG_alg.Cells[gate_idx,qbit_idx]='PHI_p45' then q_gates[gate_idx]:=qregister_calc_PHI_gate(qbit_idx-1,qbits_num,PI/4);
         if SG_alg.Cells[gate_idx,qbit_idx]='PHI_p60' then q_gates[gate_idx]:=qregister_calc_PHI_gate(qbit_idx-1,qbits_num,PI/3);
         if SG_alg.Cells[gate_idx,qbit_idx]='PHI_p90' then q_gates[gate_idx]:=qregister_calc_PHI_gate(qbit_idx-1,qbits_num,PI/2);
         if SG_alg.Cells[gate_idx,qbit_idx]='PHI_180' then q_gates[gate_idx]:=qregister_calc_PHI_gate(qbit_idx-1,qbits_num,PI);
         if SG_alg.Cells[gate_idx,qbit_idx]='PHI_m30' then q_gates[gate_idx]:=qregister_calc_PHI_gate(qbit_idx-1,qbits_num,-PI/6);
         if SG_alg.Cells[gate_idx,qbit_idx]='PHI_m45' then q_gates[gate_idx]:=qregister_calc_PHI_gate(qbit_idx-1,qbits_num,-PI/4);
         if SG_alg.Cells[gate_idx,qbit_idx]='PHI_m60' then q_gates[gate_idx]:=qregister_calc_PHI_gate(qbit_idx-1,qbits_num,-PI/3);
         if SG_alg.Cells[gate_idx,qbit_idx]='PHI_m90' then q_gates[gate_idx]:=qregister_calc_PHI_gate(qbit_idx-1,qbits_num,-PI/2);
         if SG_alg.Cells[gate_idx,qbit_idx]='SWAP1' then swap1_idx:=qbit_idx-1;
         if SG_alg.Cells[gate_idx,qbit_idx]='SWAP2' then swap2_idx:=qbit_idx-1;
         if SG_alg.Cells[gate_idx,qbit_idx]='CNOT_C' then cnot_c_idx:=qbit_idx-1;
         if SG_alg.Cells[gate_idx,qbit_idx]='CNOT_U' then cnot_u_idx:=qbit_idx-1;
      end;
      if (swap1_idx<>-1)and(swap2_idx<>-1) then q_gates[gate_idx]:=qregister_calc_SWAP_gate(swap1_idx,swap2_idx,qbits_num);
      if (cnot_c_idx<>-1)and(cnot_u_idx<>-1) then q_gates[gate_idx]:=qregister_calc_CNOT_gate(cnot_c_idx,cnot_u_idx,qbits_num);
      if length(q_gates[gate_idx])=0 then q_gates[gate_idx]:=qregister_calc_I_gate(qbits_num);
  end;

  for gate_idx:=2 to qgates_num do
  begin
       c_vector_copy(q_regs[gate_idx-1],q_regs[gate_idx]);
       qregister_apply_gate(q_regs[gate_idx],q_gates[gate_idx]);
  end;

  for gate_idx:=1 to qgates_num do
      for state_idx:=0 to states_num-1 do
      begin
          if abs(q_regs[gate_idx,state_idx].re)<1e-6 then q_regs[gate_idx,state_idx].re:=0;
          if abs(q_regs[gate_idx,state_idx].im)<1e-6 then q_regs[gate_idx,state_idx].im:=0;
      end;

  Show_all_states;
  Show_step(1);
  setlength(bin_reg,0);
end;

procedure TForm1.CB_gatesChange(Sender: TObject);
var qbit_idx:integer; gate_name:string;
begin
  CB_gates.visible:=false;
  gate_name:=CB_gates.items[CB_gates.ItemIndex];

  if (SG_alg.col=1)and((gate_name='SET')or(gate_name='RESET'))
     then SG_alg.Cells[SG_alg.col,SG_alg.row]:=gate_name;

  if (SG_alg.col>1)and
     (gate_name<>'SET')and(gate_name<>'RESET')and
     (gate_name<>'SWAP1')and(gate_name<>'SWAP2')and
     (gate_name<>'CNOT_C')and(gate_name<>'CNOT_U') then
     begin
       for qbit_idx:=1 to qbits_num do SG_alg.Cells[SG_alg.col,qbit_idx]:='';
       SG_alg.Cells[SG_alg.col,SG_alg.row]:=gate_name;
     end;

  if (SG_alg.col>1)and(gate_name='SWAP1')then
     begin
       for qbit_idx:=1 to qbits_num do
         if SG_alg.Cells[SG_alg.col,qbit_idx]<>'SWAP2' then SG_alg.Cells[SG_alg.col,qbit_idx]:='';
       SG_alg.Cells[SG_alg.col,SG_alg.row]:=gate_name;
     end;

  if (SG_alg.col>1)and(gate_name='SWAP2')then
     begin
       for qbit_idx:=1 to qbits_num do
         if SG_alg.Cells[SG_alg.col,qbit_idx]<>'SWAP1' then SG_alg.Cells[SG_alg.col,qbit_idx]:='';
       SG_alg.Cells[SG_alg.col,SG_alg.row]:=gate_name;
     end;

  if (SG_alg.col>1)and(gate_name='CNOT_C')then
     begin
       for qbit_idx:=1 to qbits_num do
         if SG_alg.Cells[SG_alg.col,qbit_idx]<>'CNOT_U' then SG_alg.Cells[SG_alg.col,qbit_idx]:='';
       SG_alg.Cells[SG_alg.col,SG_alg.row]:=gate_name;
     end;

  if (SG_alg.col>1)and(gate_name='CNOT_U')then
     begin
       for qbit_idx:=1 to qbits_num do
         if SG_alg.Cells[SG_alg.col,qbit_idx]<>'CNOT_C' then SG_alg.Cells[SG_alg.col,qbit_idx]:='';
       SG_alg.Cells[SG_alg.col,SG_alg.row]:=gate_name;
     end;
end;

procedure TForm1.FormCreate(Sender: TObject);
begin
  Reset_emul;
end;

procedure TForm1.SG_algSelectCell(Sender: TObject; aCol, aRow: Integer;
  var CanSelect: Boolean);
var R:TRect;
begin
     if (ACol>0)and(ARow>0) then
     begin
        R:=SG_alg.CellRect(ACol,ARow);
        R.Left:=R.Left+SG_alg.Left;
        R.Right:=R.Right+SG_alg.Left;
        R.top:=R.top+SG_alg.Top;
        R.Bottom:=R.Bottom+SG_alg.Top;
        CB_gates.Left:=R.Left+1;
        CB_gates.Top:=R.Top+1;
        CB_gates.width:=R.right-R.Left+1;
        CB_gates.Height:=R.Bottom-R.Top+1;
        CB_gates.Visible:=true;
     end;
     CanSelect:=true;
end;

procedure TForm1.SG_all_statesSelectCell(Sender: TObject; aCol, aRow: Integer;
  var CanSelect: Boolean);
begin
  if ACol>0 then Show_step(ACol);
end;

end.

