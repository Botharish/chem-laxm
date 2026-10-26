import './globals.css';
import './overrides.css';
import type { Metadata } from 'next';
export const metadata: Metadata = { title: 'Understanding Chemistry — Interactive Chemistry Explorer', description: 'Explore compounds, reactions, calculators and molecular models.' };
export default function RootLayout({ children }: Readonly<{ children: React.ReactNode }>) { return <html lang="en"><body>{children}</body></html>; }
