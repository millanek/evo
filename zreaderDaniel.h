//
//  zreaderDaniel.h
//  process_vcf
//
//  Created by Milan Malinsky on 28/02/2014.
//  Copyright (c) 2014 Milan Malinsky. All rights reserved.
//

#ifndef process_vcf_zreaderDaniel_h
#define process_vcf_zreaderDaniel_h

#include <zlib.h>

#include <cassert>

#include <streambuf>
#include <istream>
#include <ostream>
#include <fstream>

#include <array>

namespace dpj {
    
    /// std::streambuf wrapper around zlib, deflate and gzip streams.
    //
    class zstreambuf : public std::streambuf
    {
        static const std::size_t BUF_SIZE = 256 * 1024;
    public:
        enum class action
        {
            comp_zlib,
            decomp_zlib,
            comp_gzip,
            decomp_gzip
        };
        
    public:
        
        zstreambuf(std::streambuf* io_sbuf, action a) : mode(a), io_sbuf(io_sbuf) { init(); }
        
        void write_gzip_header(std::string const& fname, int os = 3, bool tag = false);
        
    protected:
        
        // underflow()
        //
        //   - decompresses data, fills get area.
        //
        int_type underflow()
        {
            if (gptr() < egptr())
                return traits_type::to_int_type(*gptr());
            
            std::streamsize n = inflate_buffer();
            
            setg(z_out_b, z_out_b, z_out_b + n);
            
            if (n != 0)
                return traits_type::to_int_type(*gptr());
            else
                return traits_type::eof();
        }
        
        /// overflow()
        //
        //   - compresses data, creates space in put area.
        //
        int_type overflow(int_type ch)
        {
            if (!traits_type::eq_int_type(ch, traits_type::eof()))
            {
                if (pptr() < epptr())
                    return ch;
                else
                {
                    deflate_buffer(ch);
                    setp(z_in_b, z_in_e);
                    return sputc(ch);
                }
            }
            else
            {
                setp(z_in_b, z_in_b);
                return traits_type::eof();
            }
        }
        
        /// sync() - called by std::streambuf::pubsync()
        //
        int sync()
        {
            if (mode == action::comp_zlib || mode == action::comp_gzip)
            {
                deflate_buffer(traits_type::eof());
                z_state = deflateReset(&zs);
            }
            
            if (z_state != Z_OK)
                return -1;
            
            return 0;
        }
        
    private:
        action mode;
        
        struct gz_header_info
        {
            std::string file_name, comment;
            gz_header h;
        };
        gz_header_info hdr;
        z_stream zs;
        int z_state = Z_OK;
        
        std::streambuf* io_sbuf;
        
        std::array<unsigned char, BUF_SIZE> z_in_buf;
        std::array<unsigned char, BUF_SIZE> z_out_buf;
        
        char* z_in_b, * z_in_e, * z_out_b, * z_out_e;
        
        // ///////////////////////////////////////////////////////////////////////////////
        /// zlib initialisateion routines.
        //
        void init()
        {
            z_in_b  = reinterpret_cast<char*>(z_in_buf.begin());
            z_in_e  = reinterpret_cast<char*>(z_in_buf.end());
            z_out_b = reinterpret_cast<char*>(z_out_buf.begin());
            z_out_e = reinterpret_cast<char*>(z_out_buf.end());
            
            switch (mode)
            {
                case action::decomp_zlib:
                case action::decomp_gzip:
                    setg(z_out_b, z_out_b, z_out_b);
                    z_inflate_init();
                    break;
                case action::comp_zlib:
                case action::comp_gzip:
                    setp(z_in_b, z_in_e);
                    z_deflate_init();
                    break;
                default:
                    break;
            }
        }
        
        void common_init()
        {
            zs.zalloc = Z_NULL;
            zs.zfree = Z_NULL;
            zs.opaque = Z_NULL;
            
            zs.next_in = z_in_buf.begin();
            zs.avail_in = 0;
            
            zs.next_out = z_out_buf.begin();
            zs.avail_out = static_cast<unsigned>(z_out_buf.size());
        }
        
        void z_inflate_init()
        {
            common_init();
            //
            // 15 - for the history buffer size: 2^15 == 32K, this needs to correspond to
            //      the size of the histroy buffer used in the deflate compression. Setting it too
            //      small will result in the back pointers pointing too far back. A likely
            //      Z_DATA_ERROR will result with the error message, "invalid distance too far back".
            //
            // 32 - to enable zlib and gzip decoding with automatic header detection. ***Relies on*** the
            //      header test not being done until the first call to inflate().
            //
            int r = inflateInit2(&zs, 15 + 32);
            if (r != Z_OK)
                throw std::runtime_error("zlib_streambuf::z_inflate_init() " + std::string(zError(r)));
        }
        
        void z_deflate_init()
        {
            common_init();
            
            int window_bits = mode == action::comp_gzip ? 15 + 16 : 15;
            
            int r = deflateInit2(&zs, Z_DEFAULT_COMPRESSION, Z_DEFLATED, window_bits, 8,
                                 // TODO: look into Z_RLE and Z_FILTERED.
                                 //
                                 Z_DEFAULT_STRATEGY);
            if (r != Z_OK)
                throw std::runtime_error("zstreambuf::z_deflate_init() failed.");
        }
        
        
        /// Pushes from put area to zstream.
        //
        //    - we might be called by overflow, in this case ch will be equal to either,
        //
        //    i) !traits_type::eof(), so deflate(&zs, Z_FINISH)
        //    ii) traits_type::eof(), so deflate(&zs, Z_NO_FLUSH)
        //
        //    - we might be called by sync, in this case ch == traits_type::eof().
        //
        void deflate_buffer(int_type ch)
        {
            bool is_eof = traits_type::eq_int_type(ch, traits_type::eof());
            int flush = is_eof ? Z_FINISH : Z_NO_FLUSH;
            
            // Everything in the put area is available.
            //
            zs.avail_in = static_cast<unsigned>(pptr() - pbase());
            zs.next_in = z_in_buf.begin();
            
            // We call deflate() until we _do_not_ fill the output buffer.
            //
            //   - if do not fill the output buffer then deflate() _will_ have
            //     consumed all input _or_ error.
            do
            {
                zs.next_out = z_out_buf.begin();
                zs.avail_out = static_cast<unsigned>(z_out_buf.size());
                
                z_state = deflate(&zs, flush);
                io_sbuf->sputn(z_out_b, z_out_buf.size() - zs.avail_out);
            }
            while
                (zs.avail_out == 0);
            
            if (z_state != Z_OK && z_state != Z_STREAM_END)
            {
                std::string message{"zlib_streambuf::deflate_buffer: "};
                message.append(zError(z_state));
                throw std::runtime_error{message};
            }
            
            setp(z_in_b, z_in_e);
        }
        
        /// Pulls from zstream to fill get area.
        //
        //   - Runs inflate until we get some output.
        //   - Returns zero bytes if we get Z_STREAM_END.
        //
        std::streamsize inflate_buffer()
        {
            // underflow() only calls this method when gptr() == egptr().
            //
            zs.next_out = z_out_buf.begin();
            zs.avail_out = static_cast<unsigned>(z_out_buf.size());
            
            while(z_out_buf.size() - zs.avail_out == 0)
            {
                if (z_state == Z_STREAM_END)
                {
                    inflateEnd(&zs);
                    return 0;
                }
                
                if (zs.avail_in == 0)
                {
                    zs.avail_in = static_cast<unsigned>(io_sbuf->sgetn(z_in_b, z_in_buf.size()));
                    zs.next_in = z_in_buf.begin();
                }
                z_state = inflate(&zs, Z_NO_FLUSH);
                
                if (z_state != Z_OK && z_state != Z_STREAM_END)
                {
                    std::string message{"zlib_streambuf::inflate_buffer: "};
                    message.append(zError(z_state));
                    throw std::runtime_error{message};
                }
            }
            return z_out_buf.size() - zs.avail_out; // Number of bytes written to get_buf.
        }
    };
    
    inline
    void zstreambuf::write_gzip_header(std::string const& fname, int os, bool tag)
    {
        hdr.h.text = 0;         // Agnostic.
        hdr.h.time = time(nullptr);
        hdr.h.os = os;
        hdr.h.extra = Z_NULL;
        
        hdr.file_name = fname;
        hdr.file_name.append(1, '\0');
        hdr.h.name = (uint8_t*)hdr.file_name.data();
        
        hdr.comment = "written by dpj::zstreambuf";
        hdr.comment.append(1, '\0');
        
        if (tag)
            hdr.h.comment = (uint8_t*)hdr.comment.data();
        else
            hdr.h.comment = Z_NULL;
        
        hdr.h.hcrc = 0; // no CRC at the mo.
        
        int r = deflateSetHeader(&zs, &hdr.h);
        if (!r == Z_OK)
            throw std::runtime_error("zstreambuf: couldn't write gzip header");
        
    }
    
    class zifstream : public std::istream
    {
    public:
        typedef zstreambuf::action action;
        
        zifstream(action a = action::decomp_zlib) : sb{fs.rdbuf() , a}, std::istream{&sb} { }
        
        zifstream(std::string const& s, action a = action::decomp_zlib)
        : sb{fs.rdbuf() , action::decomp_zlib}, std::istream{&sb} { fs.open(s); }
        
        zstreambuf* rdbuf() const { return const_cast<zstreambuf*>(&sb); }
        
        bool is_open() { return fs.is_open(); }
        void open(std::string const& s)
        {
            fs.open(s);
            if (fs.is_open())
                this->clear();
            else
                this->setstate(ios_base::failbit);
        }
        void close()
        {
            sb.pubsync();
            fs.close();
        }
        
        ~zifstream() { close(); }
        
    private:
        zstreambuf sb;
        std::ifstream fs;
    };
    
    
    class zofstream : public std::ostream
    {
    public:
        typedef zstreambuf::action action;
        
        zofstream(action a = action::comp_gzip) : sb{fs.rdbuf(), a} { }
        
        zofstream() : sb{fs.rdbuf() , action::comp_gzip}, std::ostream{&sb} { }
        
        zofstream(std::string const& s, action a = action::comp_gzip)
        : sb{fs.rdbuf(), a} , std::ostream{&sb}
        { fs.open(s); }
        
        zstreambuf* rdbuf() const { return const_cast<zstreambuf*>(&sb); }
        
        bool is_open() { return fs.is_open(); }
        void open(std::string const& s)
        {
            fs.open(s);
            if (fs.is_open())
                this->clear();
            else
                this->setstate(ios_base::failbit);
        }
        void close()
        {
            sb.pubsync();
            fs.close();
        }
        ~zofstream() { close(); }
        
    private:
        zstreambuf sb;
        std::ofstream fs;
    };
    
    
    
    
}


#endif
