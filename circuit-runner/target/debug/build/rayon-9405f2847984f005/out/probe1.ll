; ModuleID = 'probe1.7b1735b08454e8d2-cgu.0'
source_filename = "probe1.7b1735b08454e8d2-cgu.0"
target datalayout = "e-m:w-p270:32:32-p271:32:32-p272:64:64-i64:64-f80:128-n8:16:32:64-S128"
target triple = "x86_64-pc-windows-msvc"

%"core::iter::adapters::rev::Rev<core::iter::adapters::step_by::StepBy<core::ops::range::Range<i32>>>" = type { %"core::iter::adapters::step_by::StepBy<core::ops::range::Range<i32>>" }
%"core::iter::adapters::step_by::StepBy<core::ops::range::Range<i32>>" = type { { i32, i32 }, i64, i8, [7 x i8] }

@alloc_4aead6e2018a46d0df208d5729447de7 = private unnamed_addr constant <{ [27 x i8] }> <{ [27 x i8] c"assertion failed: step != 0" }>, align 1
@alloc_6bf223842e81ef01f9fee22d984515ee = private unnamed_addr constant <{ [89 x i8] }> <{ [89 x i8] c"/rustc/07dca489ac2d933c78d3c5158e3f43beefeb02ce\\library\\core\\src\\iter\\adapters\\step_by.rs" }>, align 1
@alloc_2fba17ab69f2dd89ea951979e54fe2fd = private unnamed_addr constant <{ ptr, [16 x i8] }> <{ ptr @alloc_6bf223842e81ef01f9fee22d984515ee, [16 x i8] c"Y\00\00\00\00\00\00\00!\00\00\00\09\00\00\00" }>, align 8

; core::iter::traits::iterator::Iterator::rev
; Function Attrs: inlinehint uwtable
define void @_ZN4core4iter6traits8iterator8Iterator3rev17ha8ab25018680b773E(ptr sret(%"core::iter::adapters::rev::Rev<core::iter::adapters::step_by::StepBy<core::ops::range::Range<i32>>>") align 8 %_0, ptr align 8 %self) unnamed_addr #0 {
start:
  call void @llvm.memcpy.p0.p0.i64(ptr align 8 %_0, ptr align 8 %self, i64 24, i1 false)
  ret void
}

; core::iter::traits::iterator::Iterator::step_by
; Function Attrs: inlinehint uwtable
define void @_ZN4core4iter6traits8iterator8Iterator7step_by17hb401ea4875245022E(ptr sret(%"core::iter::adapters::step_by::StepBy<core::ops::range::Range<i32>>") align 8 %_0, i32 %self.0, i32 %self.1, i64 %step) unnamed_addr #0 {
start:
; call core::iter::adapters::step_by::StepBy<I>::new
  call void @"_ZN4core4iter8adapters7step_by15StepBy$LT$I$GT$3new17h4e08f5d11a18362eE"(ptr sret(%"core::iter::adapters::step_by::StepBy<core::ops::range::Range<i32>>") align 8 %_0, i32 %self.0, i32 %self.1, i64 %step)
  ret void
}

; core::iter::adapters::step_by::StepBy<I>::new
; Function Attrs: inlinehint uwtable
define void @"_ZN4core4iter8adapters7step_by15StepBy$LT$I$GT$3new17h4e08f5d11a18362eE"(ptr sret(%"core::iter::adapters::step_by::StepBy<core::ops::range::Range<i32>>") align 8 %_0, i32 %iter.0, i32 %iter.1, i64 %step) unnamed_addr #0 personality ptr @__CxxFrameHandler3 {
start:
  %_7 = alloca i8, align 1
  store i8 1, ptr %_7, align 1
  %0 = icmp eq i64 %step, 0
  br i1 %0, label %bb2, label %bb1

bb2:                                              ; preds = %start
; invoke core::panicking::panic
  invoke void @_ZN4core9panicking5panic17h61d0277f5e1a7407E(ptr align 1 @alloc_4aead6e2018a46d0df208d5729447de7, i64 27, ptr align 8 @alloc_2fba17ab69f2dd89ea951979e54fe2fd) #5
          to label %unreachable unwind label %funclet_bb6

bb1:                                              ; preds = %start
  store i8 0, ptr %_7, align 1
; invoke <T as core::iter::adapters::step_by::SpecRangeSetup<T>>::setup
  %1 = invoke { i32, i32 } @"_ZN76_$LT$T$u20$as$u20$core..iter..adapters..step_by..SpecRangeSetup$LT$T$GT$$GT$5setup17hb1daf336ff7d3190E"(i32 %iter.0, i32 %iter.1, i64 %step)
          to label %bb3 unwind label %funclet_bb6

bb6:                                              ; preds = %funclet_bb6
  %2 = load i8, ptr %_7, align 1, !range !2, !noundef !3
  %3 = trunc i8 %2 to i1
  br i1 %3, label %bb5, label %bb4

funclet_bb6:                                      ; preds = %bb1, %bb2
  %cleanuppad = cleanuppad within none []
  br label %bb6

unreachable:                                      ; preds = %bb2
  unreachable

bb3:                                              ; preds = %bb1
  %iter.01 = extractvalue { i32, i32 } %1, 0
  %iter.12 = extractvalue { i32, i32 } %1, 1
  %_6 = sub i64 %step, 1
  %4 = getelementptr inbounds { i32, i32 }, ptr %_0, i32 0, i32 0
  store i32 %iter.01, ptr %4, align 8
  %5 = getelementptr inbounds { i32, i32 }, ptr %_0, i32 0, i32 1
  store i32 %iter.12, ptr %5, align 4
  %6 = getelementptr inbounds %"core::iter::adapters::step_by::StepBy<core::ops::range::Range<i32>>", ptr %_0, i32 0, i32 1
  store i64 %_6, ptr %6, align 8
  %7 = getelementptr inbounds %"core::iter::adapters::step_by::StepBy<core::ops::range::Range<i32>>", ptr %_0, i32 0, i32 2
  store i8 1, ptr %7, align 8
  ret void

bb4:                                              ; preds = %bb5, %bb6
  cleanupret from %cleanuppad unwind to caller

bb5:                                              ; preds = %bb6
  br label %bb4
}

; <T as core::iter::adapters::step_by::SpecRangeSetup<T>>::setup
; Function Attrs: inlinehint uwtable
define { i32, i32 } @"_ZN76_$LT$T$u20$as$u20$core..iter..adapters..step_by..SpecRangeSetup$LT$T$GT$$GT$5setup17hb1daf336ff7d3190E"(i32 %inner.0, i32 %inner.1, i64 %_step) unnamed_addr #0 {
start:
  %0 = insertvalue { i32, i32 } poison, i32 %inner.0, 0
  %1 = insertvalue { i32, i32 } %0, i32 %inner.1, 1
  ret { i32, i32 } %1
}

; probe1::probe
; Function Attrs: uwtable
define void @_ZN6probe15probe17hbd807a8d898518cbE() unnamed_addr #1 {
start:
  %_3 = alloca { i32, i32 }, align 4
  %_2 = alloca %"core::iter::adapters::step_by::StepBy<core::ops::range::Range<i32>>", align 8
  %_1 = alloca %"core::iter::adapters::rev::Rev<core::iter::adapters::step_by::StepBy<core::ops::range::Range<i32>>>", align 8
  store i32 0, ptr %_3, align 4
  %0 = getelementptr inbounds { i32, i32 }, ptr %_3, i32 0, i32 1
  store i32 10, ptr %0, align 4
  %1 = getelementptr inbounds { i32, i32 }, ptr %_3, i32 0, i32 0
  %2 = load i32, ptr %1, align 4, !noundef !3
  %3 = getelementptr inbounds { i32, i32 }, ptr %_3, i32 0, i32 1
  %4 = load i32, ptr %3, align 4, !noundef !3
; call core::iter::traits::iterator::Iterator::step_by
  call void @_ZN4core4iter6traits8iterator8Iterator7step_by17hb401ea4875245022E(ptr sret(%"core::iter::adapters::step_by::StepBy<core::ops::range::Range<i32>>") align 8 %_2, i32 %2, i32 %4, i64 2)
; call core::iter::traits::iterator::Iterator::rev
  call void @_ZN4core4iter6traits8iterator8Iterator3rev17ha8ab25018680b773E(ptr sret(%"core::iter::adapters::rev::Rev<core::iter::adapters::step_by::StepBy<core::ops::range::Range<i32>>>") align 8 %_1, ptr align 8 %_2)
  ret void
}

; Function Attrs: nocallback nofree nounwind willreturn memory(argmem: readwrite)
declare void @llvm.memcpy.p0.p0.i64(ptr noalias nocapture writeonly, ptr noalias nocapture readonly, i64, i1 immarg) #2

declare i32 @__CxxFrameHandler3(...) unnamed_addr #3

; core::panicking::panic
; Function Attrs: cold noinline noreturn uwtable
declare void @_ZN4core9panicking5panic17h61d0277f5e1a7407E(ptr align 1, i64, ptr align 8) unnamed_addr #4

attributes #0 = { inlinehint uwtable "target-cpu"="x86-64" }
attributes #1 = { uwtable "target-cpu"="x86-64" }
attributes #2 = { nocallback nofree nounwind willreturn memory(argmem: readwrite) }
attributes #3 = { "target-cpu"="x86-64" }
attributes #4 = { cold noinline noreturn uwtable "target-cpu"="x86-64" }
attributes #5 = { noreturn }

!llvm.module.flags = !{!0}
!llvm.ident = !{!1}

!0 = !{i32 8, !"PIC Level", i32 2}
!1 = !{!"rustc version 1.76.0 (07dca489a 2024-02-04)"}
!2 = !{i8 0, i8 2}
!3 = !{}
