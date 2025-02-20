module cudaf
   ! Interface to CUDA Runtime API

   use, intrinsic :: iso_c_binding, only: c_int

   implicit none

   integer(c_int), parameter :: cudaMemAttachGlobal = 1
   integer(c_int), parameter :: cudaMemAttachHost = 2
   integer(c_int), parameter :: cudaMemAttachSingle = 4

   integer(c_int), parameter :: cudaMemcpyHostToHost = 0
   integer(c_int), parameter :: cudaMemcpyHostToDevice = 1
   integer(c_int), parameter :: cudaMemcpyDeviceToHost = 2
   integer(c_int), parameter :: cudaMemcpyDeviceToDevice = 3
   integer(c_int), parameter :: cudaMemcpyDefault = 4

   interface

      integer(c_int) function &
             cudaMallocManaged(dPtr, size, flags) &
             bind(c, name="cudaMallocManaged")
         use, intrinsic :: iso_c_binding
         type(c_ptr), intent(out) :: dPtr
         integer(c_size_t), value :: size
         integer(c_int), value :: flags
      end function cudaMallocManaged

      integer(c_int) function &
             cudaMalloc(dPtr, size) &
             bind(c, name="cudaMalloc")
        use, intrinsic :: iso_c_binding
        type(c_ptr) :: dPtr
        integer(c_size_t), value :: size
      end function cudaMalloc

      integer(c_int) function &
             cudaMemcpy(dst, src, memSize, cpyKind) &
             bind(c, name="cudaMemcpy")
         use, intrinsic :: iso_c_binding
         type(c_ptr), value :: dst
         type(c_ptr), value :: src
         integer(c_size_t), value :: memSize
         integer(c_int), value :: cpyKind
       end function cudaMemcpy

      integer(c_int) function &
             cudaFree(dPtr) &
             bind(c, name="cudaFree")
         use, intrinsic :: iso_c_binding
         type(c_ptr), value :: dPtr
      end function cudaFree

   end interface

   contains

   ! wrapper functions

   integer function &
          device_malloc_managed(nbytes, dPtr) result(stat)
      use, intrinsic :: iso_c_binding, only: c_size_t, c_ptr
      use, intrinsic :: iso_fortran_env, only: int64
      integer(int64), intent(in) :: nbytes
      type(c_ptr), intent(inout) :: dPtr
      stat = cudaMallocManaged(dPtr, int(nbytes,c_size_t),&
                               cudaMemAttachGlobal)
   end function device_malloc_managed

   !
   integer function &
          device_malloc(nbytes, dPtr) result(stat)
      use, intrinsic :: iso_c_binding, only: c_size_t, c_ptr
      use, intrinsic :: iso_fortran_env, only: int64
      integer(int64), intent(in) :: nbytes
      type(c_ptr), intent(inout) :: dPtr
      stat = cudaMalloc(dPtr, int(nbytes,c_size_t))
   end function device_malloc

   !
   integer function &
           copy_device_to_host(nbytes, dst, src) result(stat)
      use, intrinsic :: iso_c_binding, only: c_size_t, c_ptr
      use, intrinsic :: iso_fortran_env, only: int64
      integer(int64), intent(in) :: nbytes
      type(c_ptr), intent(inout) :: dst, src
      stat = cudaMemcpy(dst, src, int(nbytes,c_size_t),&
                               cudaMemcpyDeviceToHost)
   end function copy_device_to_host

   !
   integer function &
           copy_host_to_device(nbytes, dst, src) result(stat)
      use, intrinsic :: iso_c_binding, only: c_size_t, c_ptr
      use, intrinsic :: iso_fortran_env, only: int64
      integer(int64), intent(in) :: nbytes
      type(c_ptr), intent(inout) :: dst, src
      stat = cudaMemcpy(dst, src, int(nbytes,c_size_t),&
                               cudaMemcpyHostToDevice)
   end function copy_host_to_device

   !
   integer function &
           device_free(dPtr) result(stat)
      use, intrinsic :: iso_c_binding, only: c_ptr
      type(c_ptr), intent(inout) :: dPtr
      stat = cudaFree(dPtr)
   end function device_free

end module cudaf
