	.file	"monte_carlo_SSA.c"
	.text
	.section	.rodata.str1.8,"aMS",@progbits,1
	.align 8
.LC3:
	.string	"Usage: %s <number of iterations> <output file>\n"
	.align 8
.LC21:
	.string	"total time for %d iterations:\n"
	.align 8
.LC22:
	.string	"max: %.4f, min: %.4f, avg: %.4f\n\n"
	.align 8
.LC23:
	.string	"communication time for %d iterations:\n"
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC24:
	.string	"RANK\t\t"
.LC25:
	.string	"%d         "
.LC26:
	.string	"Checkpoint %d\t"
.LC28:
	.string	"%.4fms  "
.LC29:
	.string	"w"
	.section	.rodata.str1.8
	.align 8
.LC30:
	.string	"%d, %d, %.4f, %.4f, %.4f, %d, %d\n"
	.section	.rodata.str1.1
.LC31:
	.string	"%d, %d, %d\n"
	.section	.text.startup,"ax",@progbits
	.p2align 4
	.globl	main
	.type	main, @function
main:
.LFB51:
	.cfi_startproc
	endbr64
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%r15
	pushq	%r14
	pushq	%r13
	pushq	%r12
	pushq	%rbx
	subq	$792, %rsp
	.cfi_offset 15, -24
	.cfi_offset 14, -32
	.cfi_offset 13, -40
	.cfi_offset 12, -48
	.cfi_offset 3, -56
	movl	%edi, -596(%rbp)
	movq	%rsi, -608(%rbp)
	movq	%fs:40, %rcx
	movq	%rcx, -56(%rbp)
	xorl	%ecx, %ecx
	cmpl	$3, %edi
	je	.L2
	movq	(%rsi), %rdx
	movl	$1, %edi
	xorl	%eax, %eax
	leaq	.LC3(%rip), %rsi
	call	__printf_chk@PLT
	movl	$1, %eax
.L1:
	movq	-56(%rbp), %rdx
	subq	%fs:40, %rdx
	jne	.L126
	leaq	-40(%rbp), %rsp
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	.cfi_remember_state
	.cfi_def_cfa 7, 8
	ret
.L2:
	.cfi_restore_state
	movq	8(%rsi), %rdi
	movl	$10, %edx
	xorl	%esi, %esi
	call	strtol@PLT
	leaq	-608(%rbp), %rsi
	leaq	-596(%rbp), %rdi
	movq	%rax, -776(%rbp)
	movq	-608(%rbp), %rax
	movq	16(%rax), %rax
	movq	%rax, -808(%rbp)
	call	MPI_Init@PLT
	leaq	-584(%rbp), %rsi
	leaq	ompi_mpi_comm_world(%rip), %rdi
	call	MPI_Comm_size@PLT
	leaq	-580(%rbp), %rsi
	leaq	ompi_mpi_comm_world(%rip), %rdi
	call	MPI_Comm_rank@PLT
	movslq	-584(%rbp), %rdx
	movdqa	.LC4(%rip), %xmm0
	movq	%rsp, %rcx
	movl	%edx, -600(%rbp)
	salq	$5, %rdx
	movq	%rdx, %rsi
	movaps	%xmm0, -496(%rbp)
	movq	%rdx, %rax
	andq	$-4096, %rsi
	subq	%rsi, %rcx
.L4:
	cmpq	%rcx, %rsp
	je	.L5
	subq	$4096, %rsp
	orq	$0, 4088(%rsp)
	jmp	.L4
.L5:
	andl	$4095, %eax
	subq	%rax, %rsp
	testq	%rax, %rax
	je	.L6
	orq	$0, -8(%rsp,%rax)
.L6:
	movslq	-776(%rbp), %rax
	movq	%rsp, -792(%rbp)
	movq	%rsp, %rsi
	leaq	15(,%rax,4), %rcx
	movq	%rcx, %rax
	andq	$-4096, %rcx
	andq	$-16, %rax
	subq	%rcx, %rsi
.L7:
	cmpq	%rsi, %rsp
	je	.L8
	subq	$4096, %rsp
	orq	$0, 4088(%rsp)
	jmp	.L7
.L8:
	andl	$4095, %eax
	subq	%rax, %rsp
	testq	%rax, %rax
	je	.L9
	orq	$0, -8(%rsp,%rax)
.L9:
	movq	-792(%rbp), %rdi
	xorl	%esi, %esi
	movq	%rsp, %rbx
	leaq	-448(%rbp), %r14
	movq	%rbx, -800(%rbp)
	call	memset@PLT
	xorl	%eax, %eax
	movl	$8, %ecx
	movq	%r14, %rdi
	rep stosl
	xorl	%edi, %edi
	call	time@PLT
	movl	-580(%rbp), %edi
	addl	%eax, %edi
	call	srand@PLT
	call	MPI_Wtime@PLT
	movq	-776(%rbp), %rax
	movsd	%xmm0, -816(%rbp)
	testl	%eax, %eax
	jle	.L10
	subl	$1, %eax
	movq	%rbx, -824(%rbp)
	leaq	-176(%rbp), %r15
	leaq	-480(%rbp), %r12
	movl	%eax, -828(%rbp)
	leaq	4(%rbx,%rax,4), %rax
	movq	%rax, -784(%rbp)
	movq	%rbx, -768(%rbp)
.L31:
	call	MPI_Wtime@PLT
	movq	.LC6(%rip), %rax
	xorl	%r13d, %r13d
	xorl	%ebx, %ebx
	movdqa	.LC5(%rip), %xmm3
	movsd	%xmm0, -760(%rbp)
	movq	%rax, -464(%rbp)
	movl	$20, -456(%rbp)
	movq	$0x000000000, -616(%rbp)
	movaps	%xmm3, -480(%rbp)
	.p2align 4,,10
	.p2align 3
.L11:
	movq	%r15, %rsi
	movq	%r12, %rdi
	call	prop@PLT
	movsd	-176(%rbp), %xmm2
	movapd	%xmm2, %xmm3
	addsd	-168(%rbp), %xmm3
	movsd	%xmm2, -648(%rbp)
	movapd	%xmm3, %xmm4
	addsd	-160(%rbp), %xmm4
	movsd	%xmm3, -656(%rbp)
	movapd	%xmm4, %xmm5
	addsd	-152(%rbp), %xmm5
	movsd	%xmm4, -664(%rbp)
	movapd	%xmm5, %xmm6
	addsd	-144(%rbp), %xmm6
	movsd	%xmm5, -672(%rbp)
	movapd	%xmm6, %xmm7
	addsd	-136(%rbp), %xmm7
	movsd	%xmm6, -680(%rbp)
	movapd	%xmm7, %xmm8
	addsd	-128(%rbp), %xmm8
	movsd	%xmm7, -688(%rbp)
	movsd	%xmm8, -696(%rbp)
	movapd	%xmm8, %xmm9
	addsd	-120(%rbp), %xmm9
	movsd	%xmm9, -704(%rbp)
	movapd	%xmm9, %xmm10
	addsd	-112(%rbp), %xmm10
	movsd	%xmm10, -712(%rbp)
	movapd	%xmm10, %xmm11
	addsd	-104(%rbp), %xmm11
	movsd	%xmm11, -720(%rbp)
	movapd	%xmm11, %xmm12
	addsd	-96(%rbp), %xmm12
	movsd	%xmm12, -728(%rbp)
	movapd	%xmm12, %xmm13
	addsd	-88(%rbp), %xmm13
	movsd	%xmm13, -736(%rbp)
	movapd	%xmm13, %xmm14
	addsd	-80(%rbp), %xmm14
	movsd	%xmm14, -744(%rbp)
	movapd	%xmm14, %xmm1
	addsd	-72(%rbp), %xmm1
	movapd	%xmm1, %xmm7
	addsd	-64(%rbp), %xmm7
	movsd	%xmm1, -752(%rbp)
	movsd	%xmm7, -640(%rbp)
	call	rand@PLT
	pxor	%xmm0, %xmm0
	cvtsi2sdl	%eax, %xmm0
	mulsd	.LC7(%rip), %xmm0
	call	log@PLT
	xorpd	.LC8(%rip), %xmm0
	divsd	-640(%rbp), %xmm0
	movsd	%xmm0, -624(%rbp)
	call	rand@PLT
	movsd	.LC7(%rip), %xmm0
	movsd	-640(%rbp), %xmm1
	pxor	%xmm15, %xmm15
	cvtsi2sdl	%eax, %xmm15
	movsd	-648(%rbp), %xmm2
	mulsd	%xmm1, %xmm0
	mulsd	%xmm15, %xmm0
	comisd	%xmm2, %xmm0
	jbe	.L105
	movsd	-656(%rbp), %xmm3
	comisd	%xmm3, %xmm0
	jbe	.L106
	movsd	-664(%rbp), %xmm4
	comisd	%xmm4, %xmm0
	jbe	.L107
	movsd	-672(%rbp), %xmm5
	comisd	%xmm5, %xmm0
	jbe	.L108
	movsd	-680(%rbp), %xmm6
	comisd	%xmm6, %xmm0
	jbe	.L109
	movsd	-688(%rbp), %xmm7
	comisd	%xmm7, %xmm0
	jbe	.L110
	movsd	-696(%rbp), %xmm8
	comisd	%xmm8, %xmm0
	jbe	.L111
	movsd	-704(%rbp), %xmm9
	comisd	%xmm9, %xmm0
	jbe	.L112
	movsd	-712(%rbp), %xmm10
	comisd	%xmm10, %xmm0
	jbe	.L113
	movsd	-720(%rbp), %xmm11
	comisd	%xmm11, %xmm0
	jbe	.L114
	movsd	-728(%rbp), %xmm12
	comisd	%xmm12, %xmm0
	jbe	.L115
	movsd	-736(%rbp), %xmm13
	comisd	%xmm13, %xmm0
	jbe	.L116
	movsd	-744(%rbp), %xmm14
	comisd	%xmm14, %xmm0
	jbe	.L117
	xorl	%esi, %esi
	comisd	%xmm1, %xmm0
	movsd	-752(%rbp), %xmm1
	movl	$13, %eax
	seta	%sil
	addl	$14, %esi
	comisd	%xmm1, %xmm0
	cmovbe	%eax, %esi
.L14:
	movq	%r12, %rdi
	call	update_state@PLT
	movsd	-616(%rbp), %xmm6
	leaq	(%r14,%rbx,8), %rax
	addsd	-624(%rbp), %xmm6
	movsd	.LC9(%rip), %xmm5
	comisd	%xmm6, %xmm5
	movsd	%xmm6, -616(%rbp)
	jbe	.L127
	pxor	%xmm0, %xmm0
	cvtsi2sdl	-496(%rbp,%rbx,4), %xmm0
	comisd	%xmm0, %xmm6
	jbe	.L11
	movq	%rax, -640(%rbp)
	addl	$1, %r13d
	call	MPI_Wtime@PLT
	movq	-640(%rbp), %rax
	movslq	%r13d, %rbx
	movapd	%xmm0, %xmm1
	movsd	(%rax), %xmm0
	subsd	-760(%rbp), %xmm0
	addsd	%xmm1, %xmm0
	movsd	%xmm0, (%rax)
	jmp	.L11
.L127:
	movq	%rax, -616(%rbp)
	call	MPI_Wtime@PLT
	movq	-616(%rbp), %rax
	movq	-768(%rbp), %rbx
	movapd	%xmm0, %xmm1
	movsd	(%rax), %xmm0
	subsd	-760(%rbp), %xmm0
	addq	$4, %rbx
	addsd	%xmm1, %xmm0
	movsd	%xmm0, (%rax)
	movl	-480(%rbp), %eax
	movl	%eax, -4(%rbx)
	movq	%rbx, -768(%rbp)
	cmpq	-784(%rbp), %rbx
	jne	.L31
	call	MPI_Wtime@PLT
	cmpl	$2, -828(%rbp)
	movsd	%xmm0, -616(%rbp)
	jbe	.L128
	movl	-776(%rbp), %eax
	movq	-800(%rbp), %rdi
	movdqa	.LC1(%rip), %xmm3
	movdqa	.LC2(%rip), %xmm1
	shrl	$2, %eax
	movq	-824(%rbp), %rdx
	subl	$1, %eax
	salq	$4, %rax
	leaq	16(%rdi,%rax), %rax
.L34:
	movdqu	(%rdx), %xmm0
	movdqa	%xmm1, %xmm2
	addq	$16, %rdx
	pcmpgtd	%xmm0, %xmm2
	pand	%xmm2, %xmm1
	pandn	%xmm0, %xmm2
	por	%xmm2, %xmm1
	movdqa	%xmm3, %xmm2
	pcmpgtd	%xmm0, %xmm2
	pand	%xmm2, %xmm0
	pandn	%xmm3, %xmm2
	movdqa	%xmm2, %xmm3
	por	%xmm0, %xmm3
	cmpq	%rdx, %rax
	jne	.L34
	movdqa	%xmm3, %xmm2
	movq	-776(%rbp), %rdi
	psrldq	$8, %xmm2
	movdqa	%xmm2, %xmm0
	movl	%edi, %edx
	pcmpgtd	%xmm3, %xmm0
	andl	$-4, %edx
	andb	$3, %dil
	pand	%xmm0, %xmm3
	pandn	%xmm2, %xmm0
	por	%xmm3, %xmm0
	movdqa	%xmm0, %xmm3
	psrldq	$4, %xmm3
	movdqa	%xmm3, %xmm2
	pcmpgtd	%xmm0, %xmm2
	pand	%xmm2, %xmm0
	pandn	%xmm3, %xmm2
	por	%xmm2, %xmm0
	movdqa	%xmm1, %xmm2
	psrldq	$8, %xmm2
	movd	%xmm0, %ecx
	movdqa	%xmm2, %xmm0
	pcmpgtd	%xmm1, %xmm0
	pand	%xmm0, %xmm2
	pandn	%xmm1, %xmm0
	por	%xmm2, %xmm0
	movdqa	%xmm0, %xmm2
	psrldq	$4, %xmm2
	movdqa	%xmm2, %xmm1
	pcmpgtd	%xmm0, %xmm1
	pand	%xmm1, %xmm2
	pandn	%xmm0, %xmm1
	por	%xmm2, %xmm1
	movd	%xmm1, %eax
	je	.L33
.L54:
	movq	-800(%rbp), %rdi
	movslq	%edx, %rsi
	movq	-776(%rbp), %rbx
	leaq	(%rdi,%rsi,4), %rdi
	movl	(%rdi), %esi
	cmpl	%esi, %eax
	cmovl	%esi, %eax
	cmpl	%esi, %ecx
	cmovg	%esi, %ecx
	leal	1(%rdx), %esi
	cmpl	%esi, %ebx
	jle	.L33
	movl	4(%rdi), %esi
	cmpl	%esi, %eax
	cmovl	%esi, %eax
	cmpl	%esi, %ecx
	cmovg	%esi, %ecx
	addl	$2, %edx
	cmpl	%edx, %ebx
	jle	.L33
	movl	8(%rdi), %edx
	cmpl	%edx, %eax
	cmovl	%edx, %eax
	cmpl	%edx, %ecx
	cmovg	%edx, %ecx
.L33:
	movl	%ecx, -572(%rbp)
	movl	%eax, -576(%rbp)
.L53:
	leaq	ompi_mpi_int(%rip), %r12
	leaq	-568(%rbp), %rsi
	movl	$1, %edx
	leaq	-576(%rbp), %rdi
	leaq	ompi_mpi_comm_world(%rip), %r9
	movq	%r12, %rcx
	leaq	ompi_mpi_op_max(%rip), %r8
	call	MPI_Allreduce@PLT
	leaq	ompi_mpi_op_min(%rip), %r8
	movq	%r12, %rcx
	leaq	-564(%rbp), %rsi
	leaq	-572(%rbp), %rdi
	leaq	ompi_mpi_comm_world(%rip), %r9
	movl	$1, %edx
	call	MPI_Allreduce@PLT
	movl	-568(%rbp), %r8d
	pxor	%xmm0, %xmm0
	movapd	.LC12(%rip), %xmm3
	pxor	%xmm1, %xmm1
	movapd	.LC11(%rip), %xmm2
	movl	%r8d, %eax
	subl	-564(%rbp), %eax
	cmpl	$0, -776(%rbp)
	cvtsi2sdl	%eax, %xmm0
	mulsd	.LC10(%rip), %xmm0
	cvtsi2sdl	-564(%rbp), %xmm1
	unpcklpd	%xmm0, %xmm0
	unpcklpd	%xmm1, %xmm1
	mulpd	%xmm0, %xmm3
	andpd	%xmm0, %xmm2
	addpd	%xmm1, %xmm2
	cvttpd2dq	%xmm2, %xmm2
	addpd	%xmm1, %xmm3
	cvttpd2dq	%xmm3, %xmm3
	punpcklqdq	%xmm3, %xmm2
	movapd	.LC14(%rip), %xmm3
	movaps	%xmm2, -416(%rbp)
	movapd	.LC13(%rip), %xmm2
	mulpd	%xmm0, %xmm3
	mulpd	%xmm0, %xmm2
	addpd	%xmm1, %xmm3
	addpd	%xmm1, %xmm2
	cvttpd2dq	%xmm3, %xmm3
	cvttpd2dq	%xmm2, %xmm2
	punpcklqdq	%xmm3, %xmm2
	movapd	.LC16(%rip), %xmm3
	movaps	%xmm2, -400(%rbp)
	movapd	.LC15(%rip), %xmm2
	mulpd	%xmm0, %xmm3
	mulpd	%xmm0, %xmm2
	addpd	%xmm1, %xmm3
	addpd	%xmm1, %xmm2
	cvttpd2dq	%xmm3, %xmm3
	cvttpd2dq	%xmm2, %xmm2
	punpcklqdq	%xmm3, %xmm2
	movapd	.LC18(%rip), %xmm3
	movaps	%xmm2, -384(%rbp)
	movapd	.LC17(%rip), %xmm2
	mulpd	%xmm0, %xmm3
	mulpd	%xmm0, %xmm2
	addpd	%xmm1, %xmm3
	addpd	%xmm1, %xmm2
	cvttpd2dq	%xmm3, %xmm3
	cvttpd2dq	%xmm2, %xmm2
	punpcklqdq	%xmm3, %xmm2
	movaps	%xmm2, -368(%rbp)
	movapd	.LC19(%rip), %xmm2
	mulpd	%xmm0, %xmm2
	mulpd	.LC20(%rip), %xmm0
	addpd	%xmm1, %xmm2
	addpd	%xmm1, %xmm0
	cvttpd2dq	%xmm2, %xmm1
	movdqa	%xmm1, %xmm3
	cvttpd2dq	%xmm0, %xmm0
	punpcklqdq	%xmm0, %xmm3
	movaps	%xmm3, -640(%rbp)
	movaps	%xmm3, -352(%rbp)
	jle	.L129
	movq	-776(%rbp), %rax
	movq	-800(%rbp), %rsi
	leaq	-256(%rbp), %r15
	leaq	-336(%rbp), %rdi
	leaq	-412(%rbp), %rbx
	subl	$1, %eax
	leaq	4(%rsi,%rax,4), %r9
.L42:
	movl	(%rsi), %ecx
	movq	%rdi, %rdx
	movq	%rbx, %rax
	.p2align 4,,10
	.p2align 3
.L40:
	cmpl	-4(%rax), %ecx
	jl	.L39
	cmpl	(%rax), %ecx
	jge	.L39
	addl	$1, (%rdx)
.L39:
	addq	$4, %rdx
	addq	$4, %rax
	cmpq	%r15, %rdx
	jne	.L40
	cmpl	%ecx, %r8d
	jne	.L41
	addl	$1, -260(%rbp)
.L41:
	addq	$4, %rsi
	cmpq	%r9, %rsi
	jne	.L42
.L43:
	leaq	ompi_mpi_comm_world(%rip), %rbx
	pushq	%rsi
	xorl	%r9d, %r9d
	movq	%r12, %rcx
	pushq	%rbx
	leaq	ompi_mpi_op_sum(%rip), %r8
	movl	$20, %edx
	movq	%r15, %rsi
	leaq	ompi_mpi_double(%rip), %r12
	leaq	-560(%rbp), %r13
	call	MPI_Reduce@PLT
	call	MPI_Wtime@PLT
	xorl	%r9d, %r9d
	movq	%r13, %rdi
	movq	%r12, %rcx
	movapd	%xmm0, %xmm1
	leaq	-544(%rbp), %rsi
	movq	%rbx, (%rsp)
	subsd	-816(%rbp), %xmm1
	subsd	-616(%rbp), %xmm0
	leaq	ompi_mpi_op_sum(%rip), %r8
	movl	$1, %edx
	movsd	%xmm1, -560(%rbp)
	movsd	%xmm0, -552(%rbp)
	call	MPI_Reduce@PLT
	xorl	%r9d, %r9d
	movq	%r13, %rdi
	movq	%r12, %rcx
	leaq	-536(%rbp), %rsi
	leaq	ompi_mpi_op_max(%rip), %r8
	movl	$1, %edx
	movq	%rbx, (%rsp)
	call	MPI_Reduce@PLT
	xorl	%r9d, %r9d
	movq	%r13, %rdi
	movq	%r12, %rcx
	leaq	-552(%rbp), %r13
	leaq	-528(%rbp), %rsi
	movl	$1, %edx
	movq	%rbx, (%rsp)
	leaq	ompi_mpi_op_min(%rip), %r8
	call	MPI_Reduce@PLT
	xorl	%r9d, %r9d
	movq	%r12, %rcx
	movl	$1, %edx
	leaq	-520(%rbp), %rsi
	leaq	ompi_mpi_op_sum(%rip), %r8
	movq	%r13, %rdi
	movq	%rbx, (%rsp)
	call	MPI_Reduce@PLT
	xorl	%r9d, %r9d
	movq	%r12, %rcx
	movl	$1, %edx
	leaq	-512(%rbp), %rsi
	leaq	ompi_mpi_op_max(%rip), %r8
	movq	%r13, %rdi
	movq	%rbx, (%rsp)
	call	MPI_Reduce@PLT
	xorl	%r9d, %r9d
	movq	%r12, %rcx
	movl	$1, %edx
	leaq	-504(%rbp), %rsi
	leaq	ompi_mpi_op_min(%rip), %r8
	movq	%r13, %rdi
	movq	%rbx, (%rsp)
	call	MPI_Reduce@PLT
	pushq	%rbx
	movq	-792(%rbp), %rcx
	movq	%r12, %r9
	pushq	$0
	movl	$4, %r8d
	movq	%r12, %rdx
	movl	$4, %esi
	movq	%r14, %rdi
	call	MPI_Gather@PLT
	addq	$32, %rsp
	cmpl	$0, -580(%rbp)
	je	.L130
.L44:
	call	MPI_Finalize@PLT
	xorl	%eax, %eax
	jmp	.L1
.L105:
	xorl	%esi, %esi
	jmp	.L14
.L106:
	movl	$1, %esi
	jmp	.L14
.L107:
	movl	$2, %esi
	jmp	.L14
.L108:
	movl	$3, %esi
	jmp	.L14
.L109:
	movl	$4, %esi
	jmp	.L14
.L110:
	movl	$5, %esi
	jmp	.L14
.L111:
	movl	$6, %esi
	jmp	.L14
.L112:
	movl	$7, %esi
	jmp	.L14
.L113:
	movl	$8, %esi
	jmp	.L14
.L114:
	movl	$9, %esi
	jmp	.L14
.L115:
	movl	$10, %esi
	jmp	.L14
.L116:
	movl	$11, %esi
	jmp	.L14
.L117:
	movl	$12, %esi
	jmp	.L14
.L130:
	movl	-600(%rbp), %r13d
	movl	$1, %edi
	xorl	%eax, %eax
	xorl	%ebx, %ebx
	imull	-776(%rbp), %r13d
	leaq	.LC21(%rip), %rsi
	leaq	.LC22(%rip), %r12
	movl	%r13d, %edx
	call	__printf_chk@PLT
	pxor	%xmm0, %xmm0
	movq	%r12, %rsi
	movl	$1, %edi
	cvtsi2sdl	-584(%rbp), %xmm0
	movsd	-544(%rbp), %xmm2
	movl	$3, %eax
	movsd	-528(%rbp), %xmm1
	divsd	%xmm0, %xmm2
	movsd	-536(%rbp), %xmm0
	call	__printf_chk@PLT
	movl	%r13d, %edx
	movl	$1, %edi
	xorl	%eax, %eax
	leaq	.LC23(%rip), %rsi
	call	__printf_chk@PLT
	pxor	%xmm0, %xmm0
	movq	%r12, %rsi
	movl	$1, %edi
	cvtsi2sdl	-584(%rbp), %xmm0
	movsd	-520(%rbp), %xmm2
	movsd	-504(%rbp), %xmm1
	movl	$3, %eax
	leaq	.LC25(%rip), %r12
	divsd	%xmm0, %xmm2
	movsd	-512(%rbp), %xmm0
	call	__printf_chk@PLT
	xorl	%eax, %eax
	movl	$1, %edi
	leaq	.LC24(%rip), %rsi
	call	__printf_chk@PLT
	cmpl	$0, -584(%rbp)
	jle	.L47
.L45:
	movl	%ebx, %edx
	movq	%r12, %rsi
	movl	$1, %edi
	xorl	%eax, %eax
	call	__printf_chk@PLT
	addl	$1, %ebx
	cmpl	%ebx, -584(%rbp)
	jg	.L45
.L47:
	movl	$10, %edi
	xorl	%r12d, %r12d
	leaq	.LC28(%rip), %rbx
	call	putchar@PLT
	pxor	%xmm0, %xmm0
	movsd	.LC27(%rip), %xmm1
	leaq	-496(%rbp), %rax
	cvtsi2sdl	%r13d, %xmm0
	movq	%rax, -624(%rbp)
	movl	%r13d, -648(%rbp)
	divsd	%xmm0, %xmm1
	movsd	%xmm1, -616(%rbp)
.L46:
	movq	-624(%rbp), %rax
	leaq	.LC26(%rip), %rsi
	movl	$1, %edi
	movl	(%rax,%r12,4), %edx
	xorl	%eax, %eax
	call	__printf_chk@PLT
	cmpl	$0, -584(%rbp)
	jle	.L51
	movq	-792(%rbp), %rax
	xorl	%r14d, %r14d
	leaq	(%rax,%r12,8), %r13
	.p2align 4,,10
	.p2align 3
.L50:
	movq	%rbx, %rsi
	movl	$1, %edi
	addl	$1, %r14d
	addq	$32, %r13
	movsd	-616(%rbp), %xmm0
	mulsd	-32(%r13), %xmm0
	movl	$1, %eax
	call	__printf_chk@PLT
	cmpl	%r14d, -584(%rbp)
	jg	.L50
.L51:
	movl	$10, %edi
	addq	$1, %r12
	call	putchar@PLT
	cmpq	$4, %r12
	jne	.L46
	cmpl	$0, -580(%rbp)
	movl	-648(%rbp), %r13d
	jne	.L44
	movq	-808(%rbp), %rdi
	leaq	.LC29(%rip), %rsi
	leaq	-412(%rbp), %rbx
	xorl	%r12d, %r12d
	call	fopen@PLT
	pxor	%xmm1, %xmm1
	movl	%r13d, %ecx
	movsd	-544(%rbp), %xmm0
	cvtsi2sdl	-584(%rbp), %xmm1
	movq	%rax, %r14
	pushq	%rax
	movl	-568(%rbp), %eax
	movl	-564(%rbp), %r9d
	movl	$1, %esi
	movq	%r14, %rdi
	movsd	-528(%rbp), %xmm2
	pushq	%rax
	movl	-584(%rbp), %r8d
	leaq	.LC30(%rip), %rdx
	movl	$3, %eax
	leaq	.LC31(%rip), %r13
	divsd	%xmm1, %xmm0
	movsd	-536(%rbp), %xmm1
	call	__fprintf_chk@PLT
	leaq	-416(%rbp), %r10
	movq	%rbx, -616(%rbp)
	popq	%rdx
	movq	%r10, %rbx
	popq	%rcx
.L52:
	movl	(%rbx,%r12), %ecx
	movl	(%r15,%r12), %r9d
	movq	%r13, %rdx
	movl	$1, %esi
	movq	-616(%rbp), %rax
	movq	%r14, %rdi
	movl	(%rax,%r12), %r8d
	xorl	%eax, %eax
	addq	$4, %r12
	call	__fprintf_chk@PLT
	cmpq	$76, %r12
	jne	.L52
	movl	-180(%rbp), %r9d
	movq	%r14, %rdi
	movq	%r13, %rdx
	xorl	%eax, %eax
	movl	-568(%rbp), %r8d
	movl	$1, %esi
	pshufd	$255, -640(%rbp), %xmm0
	movd	%xmm0, %ecx
	call	__fprintf_chk@PLT
	movq	%r14, %rdi
	call	fclose@PLT
	jmp	.L44
.L128:
	movl	$2147483647, %ecx
	movl	$-2147483648, %eax
	xorl	%edx, %edx
	jmp	.L54
.L129:
	leaq	-256(%rbp), %r15
	leaq	-336(%rbp), %rdi
	jmp	.L43
.L10:
	call	MPI_Wtime@PLT
	movl	$-2147483648, -576(%rbp)
	movl	$2147483647, -572(%rbp)
	movsd	%xmm0, -616(%rbp)
	jmp	.L53
.L126:
	call	__stack_chk_fail@PLT
	.cfi_endproc
.LFE51:
	.size	main, .-main
	.section	.rodata.cst16,"aM",@progbits,16
	.align 16
.LC1:
	.long	2147483647
	.long	2147483647
	.long	2147483647
	.long	2147483647
	.align 16
.LC2:
	.long	-2147483648
	.long	-2147483648
	.long	-2147483648
	.long	-2147483648
	.align 16
.LC4:
	.long	25
	.long	50
	.long	75
	.long	100
	.align 16
.LC5:
	.long	900
	.long	900
	.long	30
	.long	330
	.section	.rodata.cst8,"aM",@progbits,8
	.align 8
.LC6:
	.long	50
	.long	270
	.align 8
.LC7:
	.long	2097152
	.long	1040187392
	.section	.rodata.cst16
	.align 16
.LC8:
	.long	0
	.long	-2147483648
	.long	0
	.long	0
	.section	.rodata.cst8
	.align 8
.LC9:
	.long	0
	.long	1079574528
	.align 8
.LC10:
	.long	-1717986918
	.long	1068079513
	.section	.rodata.cst16
	.align 16
.LC11:
	.long	0
	.long	0
	.long	-1
	.long	-1
	.align 16
.LC12:
	.long	0
	.long	1073741824
	.long	0
	.long	1074266112
	.align 16
.LC13:
	.long	0
	.long	1074790400
	.long	0
	.long	1075052544
	.align 16
.LC14:
	.long	0
	.long	1075314688
	.long	0
	.long	1075576832
	.align 16
.LC15:
	.long	0
	.long	1075838976
	.long	0
	.long	1075970048
	.align 16
.LC16:
	.long	0
	.long	1076101120
	.long	0
	.long	1076232192
	.align 16
.LC17:
	.long	0
	.long	1076363264
	.long	0
	.long	1076494336
	.align 16
.LC18:
	.long	0
	.long	1076625408
	.long	0
	.long	1076756480
	.align 16
.LC19:
	.long	0
	.long	1076887552
	.long	0
	.long	1076953088
	.align 16
.LC20:
	.long	0
	.long	1077018624
	.long	0
	.long	1077084160
	.section	.rodata.cst8
	.align 8
.LC27:
	.long	0
	.long	1083129856
	.ident	"GCC: (Ubuntu 11.3.0-1ubuntu1~22.04) 11.3.0"
	.section	.note.GNU-stack,"",@progbits
	.section	.note.gnu.property,"a"
	.align 8
	.long	1f - 0f
	.long	4f - 1f
	.long	5
0:
	.string	"GNU"
1:
	.align 8
	.long	0xc0000002
	.long	3f - 2f
2:
	.long	0x3
3:
	.align 8
4:
